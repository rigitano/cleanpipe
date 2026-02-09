#!/usr/bin/env bash
set -euo pipefail

# this script calculates the total potential energy of a system. this is great for dihedral scans 
# make sure the box is big so pbc dont affect the results
#
# Usage:
#   ./energy.sh -g conf.gro -p topol.top -f charmm36
#
# Output:
#   Creates folder: energy_of_<gro_stem>/
#   Runs grompp + mdrun -rerun inside it
#   Writes RESULT.txt containing ONE NUMBER: total Potential energy (kJ/mol)

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
  potential_energy_of_<gro_stem>/RESULT.txt
EOF
}

# ---- parse args ----
GRO=""
TOP=""
FF=""
GMX_CMD="${GMX:-gmx}"
MAXWARN=1

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

# ---- output folder ----
gro_base="$(basename "$GRO")"
gro_stem="${gro_base%.*}"
outdir="potential_energy_of_${gro_stem}"
mkdir -p "$outdir"

# copy inputs so everything stays inside outdir
cp -f "$GRO" "$outdir/"
cp -f "$TOP" "$outdir/"

gro_local="$outdir/$(basename "$GRO")"
top_local="$outdir/$(basename "$TOP")"

# ---- generate MDP inside folder (single template + conditional parameters) ----
mdp_file="$outdir/energy.mdp"
cat <<EOT > "$mdp_file"
; Energy evaluation (single frame) for $FF
integrator  = md
nsteps      = 0


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

# NOTE:
# - For Martini, constraints=none and rvdw-switch=0 are effectively ignored/not used,
#   but included to keep a single MDP template with "options()" choices.
# - For "vacuum" avoid self-interaction: ensure your .gro box is large enough.

# ---- run inside outdir ----
pushd "$outdir" >/dev/null

# 1) grompp -> ener.tpr
"$GMX_CMD" grompp \
  -f "$(basename "$mdp_file")" \
  -c "$(basename "$gro_local")" \
  -p "$(basename "$top_local")" \
  -o ener.tpr \
  -maxwarn "$MAXWARN"

# 2) mdrun -rerun -> ener.edr
"$GMX_CMD" mdrun -s ener.tpr -rerun "$(basename "$gro_local")" -deffnm ener

# 3) extract Potential term to potential.xvg
# We select by name; this usually works. If your build requires a numeric index,
# run once interactively:  gmx energy -f ener.edr  (see the index for Potential)
echo "Potential" | "$GMX_CMD" energy -f ener.edr -o potential.xvg >/dev/null 2>&1

# 4) write RESULT.txt with a single number (kJ/mol)
# Grab 2nd column from first non-comment line
pot="$(awk '!/^[@#]/{print $2; exit}' potential.xvg)"

printf "%s\n" "$pot" > RESULT.txt

popd >/dev/null

echo "Done."
echo "Folder: $outdir"
echo "Potential energy (kJ/mol): $(cat "$outdir/RESULT.txt")"
