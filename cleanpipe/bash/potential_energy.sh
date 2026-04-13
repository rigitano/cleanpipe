#!/usr/bin/env bash
set -euo pipefail

# This script calculates the single-point POTENTIAL energy of a structure
# without minimization and without changing coordinates.
#
# It writes:
#   dihedral_scan_gromacs/potential_energy_of_<gro_stem>/
# containing:
#   energy.mdp
#   ener.tpr
#   ener.edr
#   ener.log
#   potential.xvg
#   RESULT.txt
#
# RESULT.txt contains ONE NUMBER in kJ/mol.
# The same number is also printed to stdout.

usage() {
  cat <<'EOF'
Usage:
  energy.sh -g <conf.gro> -p <topol.top> -f <martini3|charmm36> [-x <gmx>] [-w <maxwarn>]

Options:
  -g   input coordinate file (.gro recommended)
  -p   input topology (.top)
  -f   force field family: martini3 or charmm36
  -x   gmx command (default: $GMX if set, else "gmx")
  -w   maxwarn for grompp (default: 2)
  -h   show this help

Writes:
  dihedral_scan_gromacs/potential_energy_of_<gro_stem>/RESULT.txt
EOF
}

GRO=""
TOP=""
FF=""
GMX_CMD="${GMX:-gmx}"
MAXWARN=2

while getopts ":g:p:f:x:w:h" opt; do
  case "$opt" in
    g) GRO="$OPTARG" ;;
    p) TOP="$OPTARG" ;;
    f) FF="$OPTARG" ;;
    x) GMX_CMD="$OPTARG" ;;
    w) MAXWARN="$OPTARG" ;;
    h) usage; exit 0 ;;
    \?) echo "ERROR: Unknown option: -$OPTARG" >&2; usage; exit 2 ;;
    :)  echo "ERROR: Missing argument for -$OPTARG" >&2; usage; exit 2 ;;
  esac
done

if [[ -z "$GRO" || -z "$TOP" || -z "$FF" ]]; then
  echo "ERROR: -g, -p, and -f are required." >&2
  usage
  exit 2
fi

if [[ "$FF" != "martini3" && "$FF" != "charmm36" ]]; then
  echo "ERROR: -f must be 'martini3' or 'charmm36'." >&2
  exit 2
fi

if [[ ! -f "$GRO" ]]; then
  echo "ERROR: Coordinate file not found: $GRO" >&2
  exit 2
fi

if [[ ! -f "$TOP" ]]; then
  echo "ERROR: Topology file not found: $TOP" >&2
  exit 2
fi

# Absolute paths keep topology #includes working when running inside outdir
GRO_ABS="$(readlink -f "$GRO")"
TOP_ABS="$(readlink -f "$TOP")"

options() {
  declare -A map
  for pair in "$@"; do
    key="${pair%%=*}"
    val="${pair#*=}"
    map["$key"]="$val"
  done
  printf '%s\n' "${map[$FF]}"
}

scan_root="dihedral_scan_gromacs2"
mkdir -p "$scan_root"

gro_base="$(basename "$GRO_ABS")"
gro_stem="${gro_base%.*}"
outdir="${scan_root}/potential_energy_of_${gro_stem}"
mkdir -p "$outdir"

# Optional symlinks for convenience
ln -sfn "$GRO_ABS" "$outdir/$(basename "$GRO_ABS")"
ln -sfn "$TOP_ABS" "$outdir/$(basename "$TOP_ABS")"

mdp_file="$outdir/energy.mdp"
cat > "$mdp_file" <<EOT
; Single-point potential energy evaluation for $FF
integrator               = md
nsteps                   = 0
dt                       = $(options charmm36=0.002 martini3=0.02)

; Neighbor searching
cutoff-scheme            = Verlet
nstlist                  = $(options charmm36=10 martini3=20)
rlist                    = $(options charmm36=1.2 martini3=1.1)
pbc                      = xyz

; Electrostatics
coulombtype              = $(options charmm36=PME martini3=Reaction-Field)
rcoulomb                 = $(options charmm36=1.2 martini3=1.1)
epsilon_r                = $(options charmm36=1 martini3=15)
epsilon_rf               = $(options charmm36=0 martini3=0)

; van der Waals
vdwtype                  = Cut-off
rvdw                     = $(options charmm36=1.2 martini3=1.1)
vdw-modifier             = $(options charmm36=Force-switch martini3=Potential-shift-verlet)
rvdw-switch              = $(options charmm36=1.0 martini3=0)
DispCorr                 = $(options charmm36=no martini3=no)

; Constraints
constraints              = $(options charmm36=h-bonds martini3=none)

; No thermostat/barostat
tcoupl                   = no
pcoupl                   = no

; Output
nstenergy                = 1
nstlog                   = 1
nstxout                  = 0
nstvout                  = 0
nstfout                  = 0
nstxout-compressed       = 0
EOT

pushd "$outdir" >/dev/null

"$GMX_CMD" grompp \
  -f "$(basename "$mdp_file")" \
  -c "$GRO_ABS" \
  -p "$TOP_ABS" \
  -o ener.tpr \
  -maxwarn "$MAXWARN"

"$GMX_CMD" mdrun \
  -s ener.tpr \
  -rerun "$GRO_ABS" \
  -deffnm ener

printf "Potential\n0\n" | "$GMX_CMD" energy -f ener.edr -o potential.xvg >/dev/null 2>&1

# Use the LAST data line, not the first one
pot="$(awk '!/^[@#]/ {val=$2} END {if (val=="") exit 1; print val}' potential.xvg)"

printf "%s\n" "$pot" > RESULT.txt

popd >/dev/null

printf "%s\n" "$pot"