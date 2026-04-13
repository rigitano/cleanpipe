#!/usr/bin/env bash
set -euo pipefail

# This script calculates the total potential energy of a system (great for dihedral scans).
# It writes a folder:
#   dihedral_scan_gromacs/potential_energy_of_<gro_stem>/
# containing RESULT.txt with ONE NUMBER (kJ/mol).

usage() {
  cat <<'EOF'
Usage:
  energy.sh -g <conf.gro> -p <topol.top> -f <martini3|charmm36> [-x <gmx>] [-w <maxwarn>]

Options:
  -g   input .gro (or .pdb, but .gro recommended)
  -p   input .top
  -f   forcefield: martini3 or charmm36
  -x   gmx command (default: $GMX if set, else "gmx")
  -w   maxwarn for grompp (default: 1)

Writes:
  dihedral_scan_gromacs/potential_energy_of_<gro_stem>/RESULT.txt
EOF
}

# ---- parse args ----
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
    \?) echo "Unknown option: -$OPTARG" >&2; usage; exit 2 ;;
    :)  echo "Missing argument for -$OPTARG" >&2; usage; exit 2 ;;
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
  echo "ERROR: GRO file not found: $GRO" >&2
  exit 2
fi
if [[ ! -f "$TOP" ]]; then
  echo "ERROR: TOP file not found: $TOP" >&2
  exit 2
fi

# ---- absolute paths so running in outdir doesn't break relative #includes ----
GRO_ABS="$(readlink -f "$GRO")"
TOP_ABS="$(readlink -f "$TOP")"
TOP_DIR="$(dirname "$TOP_ABS")"

# ---- your requested "options()" pattern ----
options() {
  declare -A map
  for pair in "$@"; do
    key="${pair%%=*}"
    val="${pair#*=}"
    map["$key"]="$val"
  done
  echo "${map[$FF]}"
}

# ---- parent folder for all scan points ----
scan_root="dihedral_scan_gromacs"
mkdir -p "$scan_root"

# ---- output folder (inside dihedral_scan_gromacs/) ----
gro_base="$(basename "$GRO_ABS")"
gro_stem="${gro_base%.*}"
outdir="${scan_root}/potential_energy_of_${gro_stem}"
mkdir -p "$outdir"

# (optional) keep references to inputs in the outdir without copying (copying can break #includes)
ln -sf "$GRO_ABS" "$outdir/$(basename "$GRO_ABS")"
ln -sf "$TOP_ABS" "$outdir/$(basename "$TOP_ABS")"

# ---- generate MDP inside folder (single template + conditional parameters) ----
mdp_file="$outdir/energy.mdp"
cat <<EOT > "$mdp_file"
; Energy evaluation (single frame) for $FF
integrator  = md
nsteps      = 0

; in case you want to minimize first
; Integrator =	steep
; emtol      =	1000 ;100
; emstep     =	0.01
; nsteps     =    100000 ; this is the max value to be used just if emtol is never reached

; Neighbor searching
cutoff-scheme = Verlet
nstlist     = $(options charmm36=10 martini3=20)
rlist       = $(options charmm36=1.2 martini3=1.1)

; Electrostatics
rcoulomb    = $(options charmm36=1.2 martini3=1.1)
coulombtype = $(options charmm36=PME martini3=reaction-field)

epsilon_r   = $(options charmm36=1 martini3=15)
epsilon_rf  = 0

; Van der Waals
rvdw        = $(options charmm36=1.2 martini3=1.1)
vdwtype     = cutoff
vdw-modifier= $(options charmm36=force-switch martini3=Potential-shift-verlet)
rvdw-switch = $(options charmm36=1.0 martini3=0)
DispCorr    = no

constraints = $(options charmm36=h-bonds martini3=none)

tcoupl      = no
pcoupl      = no
EOT

# ---- run inside outdir ----
pushd "$outdir" >/dev/null

# 1) grompp -> ener.tpr
"$GMX_CMD" grompp \
  -f "$(basename "$mdp_file")" \
  -c "$GRO_ABS" \
  -p "$TOP_ABS" \
  -o ener.tpr \
  -maxwarn "$MAXWARN"


# 2) mdrun -rerun -> ener.edr
"$GMX_CMD" mdrun -s ener.tpr -rerun "$GRO_ABS" -deffnm ener

# 3) extract Potential term to potential.xvg
echo "Potential" | "$GMX_CMD" energy -f ener.edr -o potential.xvg >/dev/null 2>&1

# 4) write RESULT.txt with a single number (kJ/mol)
pot="$(awk '!/^[@#]/{print $2; exit}' potential.xvg)"
printf "%s\n" "$pot" > RESULT.txt

popd >/dev/null

echo "Done."
echo "Folder: $outdir"
echo "Potential energy (kJ/mol): $(cat "$outdir/RESULT.txt")"
