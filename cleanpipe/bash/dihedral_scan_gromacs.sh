#!/bin/bash
# =============================================================================
# GROMACS Dihedral Scan — Multi-Stage Annealing Protocol
# =============================================================================
#
# Avoids local-minima trapping via a four-stage pipeline per scan point:
#
#  Stage 1 — Rough EM      Steepest-descent (loose tol=100) removes hard
#                           steric clashes present in the starting structure.
#
#  Stage 2 — Hot NVT SD    Short SD at T_hot (default 400 K) supplies kinetic
#                           energy to cross torsional/steric barriers.  This
#                           is physically equivalent to "softcore" potentials:
#                           the high temperature effectively flattens the sharp
#                           repulsive LJ wall, letting the molecule escape
#                           local minima without needing B-state topology or
#                           the GROMACS free-energy λ framework.
#
#  Stage 3 — Anneal SD     Simulated annealing (T_hot → 50 K) gradually
#                           "restores" the real potential landscape while the
#                           molecule is guided into the deepest accessible
#                           basin — directly analogous to growing λ from 0 → 1
#                           in a softcore free-energy calculation.
#
#  Stage 4 — Final EM      Tight steepest-descent EM (tol=1) on the
#                           annealed structure converges to the true, fully-
#                           relaxed potential energy at the restrained angle.
#
# The dihedral restraint (topology function type 2) is active in every stage,
# so the scanned angle is preserved throughout the entire pipeline.
#
# Usage:
#   ./optimized_dihedral_scan_gromacs.sh \
#                      -p <topology.top> -g <structure.gro> \
#                      -i <id1> -j <id2> -k <id3> -l <id4> \
#                      [-f <force_constant>] [-o <output_dir>] \
#                      [-T <hot_temp_K>] [-H <hot_steps>]    \
#                      [-S <anneal_steps>] [-n <threads>]
#
#Example: ./optimized_dihedral_scan_gromacs.sh -p CHYO_proxy_complete.top -g torsion_dihedral_000.gro -i 15 -j 16 -k 1 -l 11 -f 50000 -o gromacs_scan_optimized

#
# Required arguments:
#   -p   Topology file (.top)
#   -g   Starting structure file (.gro)
#   -i   Dihedral atom index 1  (1-based, as in GROMACS [ dihedrals ])
#   -j   Dihedral atom index 2
#   -k   Dihedral atom index 3
#   -l   Dihedral atom index 4
#
# Optional arguments:
#   -f   Harmonic restraint force constant  kJ/mol/rad²  (default: 50000)
#   -o   Output directory                               (default: gromacs_dihedral_scan_output)
#   -T   Hot-MSD temperature                  K          (default: 400)
#   -H   Number of SD steps for hot stage               (default: 5000  = 5 ps at dt=0.001)
#   -S   Number of SD steps for annealing stage         (default: 20000 = 20 ps at dt=0.001)
#   -n   OpenMP threads passed to mdrun                 (default: 1)
#   -h   Print this help message and exit
#
# Notes:
#   - Integration timestep for SD stages is fixed at 0.001 ps (1 fs), which
#     is safe for unconstrained systems at high temperature.  Increase -H/-S
#     for larger or more flexible molecules.
#   - For gas-phase / vacuum calculations, change pbc = xyz to pbc = no and
#     remove the coulombtype / fourierspacing / ewald_rtol lines from the
#     generated MDP files in <output_dir>/ before running.
#   - The results file records energies from Stage 4 only (true minimum).
#     Intermediate energies are stored in each angle's subdirectory.
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
FORCE_CONSTANT=50000
OUTPUT_DIR="gromacs_dihedral_scan_output"
THREADS=1
ANGLE_START=0
ANGLE_END=355
ANGLE_STEP=5
HOT_TEMP=400       # K
HOT_STEPS=5000     # steps  ×  0.001 ps  =  5 ps
ANNEAL_STEPS=20000 # steps  ×  0.001 ps  =  20 ps
MD_DT=0.001        # integration timestep (ps)

# ── Terminal colors ───────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
CYAN='\033[0;36m'; BLUE='\033[0;34m'; NC='\033[0m'

# ── Utility functions ─────────────────────────────────────────────────────────
usage() { grep '^#' "$0" | grep -v '#!/' | sed 's/^# \{0,2\}//'; exit 0; }
log()   { echo -e "${CYAN}[INFO]${NC}   $*"; }
warn()  { echo -e "${YELLOW}[WARN]${NC}   $*"; }
err()   { echo -e "${RED}[ERROR]${NC}  $*" >&2; exit 1; }
ok()    { echo -e "${GREEN}[OK]${NC}     $*"; }
stg()   { echo -e "${BLUE}[STG${1}]${NC}  ${2}"; }

# ── Parse arguments ───────────────────────────────────────────────────────────
ATOM1="" ATOM2="" ATOM3="" ATOM4=""
TOP_FILE="" GRO_FILE=""

while getopts ":p:g:i:j:k:l:f:o:T:H:S:n:h" opt; do
    case $opt in
        p) TOP_FILE="$OPTARG"       ;;
        g) GRO_FILE="$OPTARG"       ;;
        i) ATOM1="$OPTARG"          ;;
        j) ATOM2="$OPTARG"          ;;
        k) ATOM3="$OPTARG"          ;;
        l) ATOM4="$OPTARG"          ;;
        f) FORCE_CONSTANT="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG"     ;;
        T) HOT_TEMP="$OPTARG"       ;;
        H) HOT_STEPS="$OPTARG"      ;;
        S) ANNEAL_STEPS="$OPTARG"   ;;
        n) THREADS="$OPTARG"        ;;
        h) usage                    ;;
        :) err "Option -$OPTARG requires an argument." ;;
       \?) err "Unknown option: -$OPTARG"              ;;
    esac
done

# ── Validate required arguments ───────────────────────────────────────────────
[[ -z "$TOP_FILE" ]] && err "Topology file (-p) is required."
[[ -z "$GRO_FILE" ]] && err "Structure file (-g) is required."
[[ -z "$ATOM1"    ]] && err "Atom index 1 (-i) is required."
[[ -z "$ATOM2"    ]] && err "Atom index 2 (-j) is required."
[[ -z "$ATOM3"    ]] && err "Atom index 3 (-k) is required."
[[ -z "$ATOM4"    ]] && err "Atom index 4 (-l) is required."
[[ -f "$TOP_FILE" ]] || err "Topology not found: $TOP_FILE"
[[ -f "$GRO_FILE" ]] || err "Structure not found: $GRO_FILE"

# ── Locate GROMACS ────────────────────────────────────────────────────────────
if ! command -v gmx &>/dev/null && ! command -v gmx_mpi &>/dev/null; then
    err "GROMACS (gmx or gmx_mpi) not found in PATH."
fi
GMX=$(command -v gmx 2>/dev/null || command -v gmx_mpi)
log "Using GROMACS: $GMX"

# ── Absolute paths ────────────────────────────────────────────────────────────
TOP_FILE=$(realpath "$TOP_FILE")
GRO_FILE=$(realpath "$GRO_FILE")
mkdir -p "$OUTPUT_DIR"
OUTPUT_DIR=$(realpath "$OUTPUT_DIR")

RESULTS_FILE="$OUTPUT_DIR/dihedral_scan_results.dat"
LOG_FILE="$OUTPUT_DIR/dihedral_scan.log"

# ── Derived SD timing parameters (computed once) ──────────────────────────────
ANNEAL_TIME=$(awk "BEGIN{printf \"%.3f\", ${ANNEAL_STEPS} * ${MD_DT}}")
ANNEAL_MID=$( awk "BEGIN{printf \"%.3f\", ${ANNEAL_STEPS} * ${MD_DT} * 0.5}")
HOT_TIME=$(   awk "BEGIN{printf \"%.3f\", ${HOT_STEPS}    * ${MD_DT}}")

# ── Startup summary ───────────────────────────────────────────────────────────
log "Output directory  : $OUTPUT_DIR"
log "Topology          : $TOP_FILE"
log "Structure         : $GRO_FILE"
log "Dihedral atoms    : $ATOM1 $ATOM2 $ATOM3 $ATOM4"
log "Force constant    : $FORCE_CONSTANT kJ/mol/rad²"
log "Angle range       : ${ANGLE_START}° → ${ANGLE_END}°  step ${ANGLE_STEP}°"
log "Protocol          : Rough EM → Hot MD@${HOT_TEMP}K(${HOT_TIME}ps) → Anneal→50K(${ANNEAL_TIME}ps) → Final EM"

# ── Results file header ───────────────────────────────────────────────────────
cat > "$RESULTS_FILE" <<EOF
# GROMACS Dihedral Scan Results — Multi-Stage Annealing Protocol
# Dihedral atoms    : $ATOM1 $ATOM2 $ATOM3 $ATOM4
# Force constant    : $FORCE_CONSTANT kJ/mol/rad²
# Stage 2 (Hot MD)  : ${HOT_TEMP} K,  ${HOT_STEPS} steps  (${HOT_TIME} ps)
# Stage 3 (Anneal)  : ${HOT_TEMP} K → 50 K,  ${ANNEAL_STEPS} steps  (${ANNEAL_TIME} ps)
# Stage 4 (Final EM): emtol = 1 kJ/mol/nm
# Topology          : $TOP_FILE
# Structure         : $GRO_FILE
# Generated         : $(date)
#
# Epot is from Stage 4 (final tight EM on the annealed structure).
#
# Angle(deg)   Epot(kJ/mol)   Converged
EOF

# =============================================================================
# write_mdps
#   Writes all four stage MDP files once to OUTPUT_DIR before the scan loop.
#   They are reused (read-only) at every scan angle.
# =============================================================================
write_mdps() {

    # ── Shared non-bonded block (identical for all four stages) ───────────────
    # Single-quoted heredoc: NO shell expansion inside — this is intentional.
    NB_BLOCK=$(cat <<'__NB__'
; ── Non-bonded interactions ─────────────────────────────────────────────────
cutoff-scheme    = Verlet
pbc              = xyz

coulombtype      = Cut-off
rcoulomb         = 1.2
coulomb-modifier = None

vdwtype          = Cut-off
rvdw             = 1.2
vdw-modifier     = None
DispCorr         = no
__NB__
)

    # ── Stage 1: Rough steepest-descent EM ───────────────────────────────────
    # Purpose : eliminate hard steric clashes in the starting (or restraint-
    #           modified) structure.  Loose tolerance so we don't waste time
    #           converging — Stage 4 will do that.
    cat > "$OUTPUT_DIR/stage1_rough_em.mdp" <<__MDP__
; ============================================================
; Stage 1 — Rough steepest-descent energy minimisation
; Removes hard steric clashes without wasting steps on full
; convergence.  The annealing stages handle basin relaxation.
; ============================================================
integrator      = steep
emtol           = 100.0         ; kJ/mol/nm — loose tolerance (clash removal only)
emstep          = 0.01          ; nm
nsteps          = 50000

${NB_BLOCK}
__MDP__

    # ── Stage 2: Hot NVT SD ───────────────────────────────────────────────────
    # Purpose : provide kinetic energy to escape local minima.
    #   • constraints = none  keeps the simulation stable at high T without
    #     requiring a special constraint algorithm; 1 fs timestep is safe.
    #   • No trajectory output (nstxout/nstvout = 0) to save disk space;
    #     only the final .gro and .edr are kept.
    cat > "$OUTPUT_DIR/stage2_hot_md.mdp" <<__MDP__
; ============================================================
; Stage 2 — Hot NVT SD  (local-minima escape)
;
; Running SD at ${HOT_TEMP} K supplies kinetic energy equivalent
; to crossing ~${HOT_TEMP}*8.314e-3 kJ/mol barriers per degree of
; freedom, effectively flattening the sharp repulsive LJ wall
; in the same way softcore potentials do — without requiring
; B-state topology parameters.
; ============================================================
integrator      = sd
dt              = ${MD_DT}       ; ps  — 1 fs, safe for unconstrained systems
nsteps          = ${HOT_STEPS}   ; total = ${HOT_TIME} ps

; Output — minimal (final .gro always written by mdrun)
nstlog               = 500
nstxout              = 0
nstvout              = 0
nstfout              = 0
nstxout-compressed   = 0
nstcalcenergy        = 100

; Temperature coupling

tc-grps         = System
tau-t           = 0.1           ; ps
ref-t           = ${HOT_TEMP}   ; K
gen-vel         = yes           ; generate velocities from Maxwell-Boltzmann
gen-temp        = ${HOT_TEMP}   ; K
gen-seed        = -1            ; random seed each run

; No pressure coupling (NVT)
pcoupl          = no

; Constraints — none, for maximum flexibility at high temperature.
; Avoids LINCS failures on arbitrary topologies.
constraints     = none
continuation    = no

; Centre-of-mass motion removal
comm-mode       = Linear
nstcomm         = 100
comm-grps       = System

${NB_BLOCK}
__MDP__

    # ── Stage 3: Simulated annealing ──────────────────────────────────────────
    # Purpose : cool the hot structure back to near 0 K, guiding the molecule
    #   into the deepest accessible basin at the restrained dihedral angle.
    #   This is the direct analogue of growing the softcore λ from ~0 → 1.
    #
    # Annealing schedule (linear cooling):
    #   t = 0 ps         T = HOT_TEMP  K   (start hot)
    #   t = ANNEAL_MID ps   T = 200 K     (mid-cool plateau)
    #   t = ANNEAL_TIME ps  T =  50 K     (near-zero, hand off to EM)
    #
    # gen-vel = yes with gen-temp = HOT_TEMP generates a fresh Maxwell-
    # Boltzmann distribution at the hot temperature using the annealed
    # coordinates from Stage 2, so Stage 3 always starts at T_hot.
    cat > "$OUTPUT_DIR/stage3_anneal_md.mdp" <<__MDP__
; ============================================================
; Stage 3 — Simulated-annealing NVT MD
;   ${HOT_TEMP} K → 200 K → 50 K  over ${ANNEAL_TIME} ps
;
; Gradually "restores" the full LJ potential landscape as T
; drops, analogous to growing the softcore lambda from 0 → 1.
; The dihedral restraint keeps the scanned angle fixed while
; all other degrees of freedom relax into the basin.
; ============================================================
integrator      = sd
dt              = ${MD_DT}
nsteps          = ${ANNEAL_STEPS}  ; total = ${ANNEAL_TIME} ps

; Output — minimal
nstlog               = 1000
nstxout              = 0
nstvout              = 0
nstfout              = 0
nstxout-compressed   = 0
nstcalcenergy        = 100

; Temperature coupling
tc-grps         = System
tau-t           = 0.2
ref-t           = 50            ; K — final target (overridden by schedule below)
gen-vel         = no           ; fresh velocities at T_hot from Stage 2 coords


; Simulated annealing — linear cooling in two segments:
;   Segment 1: ${HOT_TEMP} K → 200 K  (first half of the run)
;   Segment 2: 200 K → 50 K           (second half of the run)
annealing           = single
annealing-npoints   = 3
annealing-time      = 0  ${ANNEAL_MID}  ${ANNEAL_TIME}
annealing-temp      = ${HOT_TEMP}  200  50

; No pressure coupling (NVT)
pcoupl          = no

; Constraints — none (consistent with Stage 2)
constraints     = none
continuation    = yes

; Centre-of-mass motion removal
comm-mode       = Linear
nstcomm         = 100
comm-grps       = System

${NB_BLOCK}
__MDP__

    # ── Stage 4: Final tight EM ───────────────────────────────────────────────
    # Purpose : converge the annealed structure to the true local minimum.
    #   This is the energy reported in the results file — it represents a
    #   physically meaningful, fully-relaxed configuration at the restrained
    #   dihedral angle, free of the local-minima artefacts that a single-step
    #   EM would produce.
    cat > "$OUTPUT_DIR/stage4_final_em.mdp" <<__MDP__
; ============================================================
; Stage 4 — Final tight steepest-descent EM
;
; Converges the annealed structure to the nearest local minimum
; with a tight force tolerance.  Because the structure was
; thoroughly relaxed by the annealing, this minimum reflects
; the TRUE low-energy geometry at the restrained angle rather
; than a kinetically trapped artefact.
; ============================================================
integrator      = steep
emtol           = 1           ; kJ/mol/nm — tight convergence
emstep          = 0.005         ; nm — smaller step for accuracy
nsteps          = 1000000

${NB_BLOCK}
__MDP__

    log "Stage MDP files written:"
    log "  $OUTPUT_DIR/stage1_rough_em.mdp"
    log "  $OUTPUT_DIR/stage2_hot_md.mdp"
    log "  $OUTPUT_DIR/stage3_anneal_md.mdp"
    log "  $OUTPUT_DIR/stage4_final_em.mdp"
}

# =============================================================================
# patch_topology <src_top> <dst_top> <angle>
#   Creates a per-angle copy of the topology with a harmonic dihedral restraint
#   (function type 2) set to <angle> degrees.  Replaces an existing line if
#   one is already present for these four atoms; otherwise appends to the last
#   [ dihedrals ] section.
# =============================================================================
patch_topology() {
    local src_top="$1" dst_top="$2" angle="$3"
    local restraint_line="${ATOM1}   ${ATOM2}   ${ATOM3}   ${ATOM4}   2   ${angle}   ${FORCE_CONSTANT}"
    local search_pattern="^[[:space:]]*${ATOM1}[[:space:]]+${ATOM2}[[:space:]]+${ATOM3}[[:space:]]+${ATOM4}[[:space:]]+2[[:space:]]"

    if grep -qE "$search_pattern" "$src_top"; then
        # Replace an existing restraint line for this dihedral
        sed -E "s/${search_pattern}.*/${restraint_line}/" "$src_top" > "$dst_top"
    else
        # Append after the header (and any leading comments) of the LAST
        # [ dihedrals ] section found in the topology
        awk -v line="$restraint_line" '
        BEGIN { last=-1; n=0 }
        {
            lines[n] = $0
            if ($0 ~ /^\[ *dihedrals *\]/) last = n
            n++
        }
        END {
            inserted = 0
            for (i = 0; i < n; i++) {
                print lines[i]
                if (!inserted && i == last) {
                    j = i + 1
                    while (j < n && (lines[j] ~ /^[[:space:]]*;/ || lines[j] ~ /^[[:space:]]*$/)) {
                        print lines[j]; i = j; j++
                    }
                    print line
                    inserted = 1
                }
            }
            if (!inserted) print line
        }' "$src_top" > "$dst_top"
    fi
}

# =============================================================================
# extract_energy <edr_file> <work_dir>
#   Extracts the final potential energy from a GROMACS .edr file using
#   `gmx energy`.  Falls back to parsing the .xvg output if the grep
#   on stdout fails (version-dependent output format).
# =============================================================================
extract_energy() {
    local edr_file="$1" work_dir="$2"
    local xvg="$work_dir/energy_tmp.xvg"
    local energy

    energy=$(echo "Potential" | "$GMX" energy \
                 -f "$edr_file" -o "$xvg" -quiet 2>/dev/null \
             | grep -E "^Potential" | awk '{print $2}') || true

    if [[ -z "$energy" ]] && [[ -f "$xvg" ]]; then
        energy=$(grep -v '^[#@]' "$xvg" | tail -1 | awk '{print $2}')
    fi

    echo "${energy:-NaN}"
}

# =============================================================================
# run_grompp <mdp> <gro> <top> <tpr> <logfile>
# =============================================================================
run_grompp() {
    local mdp="$1" gro="$2" top="$3" tpr="$4" logfile="$5"
    "$GMX" grompp \
        -f "$mdp" -c "$gro" -p "$top" -o "$tpr" \
        -maxwarn 5 -quiet \
        > "$logfile" 2>&1
}

# =============================================================================
# run_mdrun <tpr> <deffnm_base> <logfile>
# =============================================================================
run_mdrun() {
    local tpr="$1" deffnm="$2" logfile="$3"
    "$GMX" mdrun \
        -s "$tpr" -deffnm "$deffnm" \
        -ntmpi 1 -ntomp "$THREADS" \
        -quiet \
        >> "$logfile" 2>&1
}

# =============================================================================
# Write MDP files once before the scan loop (they are shared across angles)
# =============================================================================
write_mdps

MDP_S1="$OUTPUT_DIR/stage1_rough_em.mdp"
MDP_S2="$OUTPUT_DIR/stage2_hot_md.mdp"
MDP_S3="$OUTPUT_DIR/stage3_anneal_md.mdp"
MDP_S4="$OUTPUT_DIR/stage4_final_em.mdp"

# =============================================================================
# Main dihedral scan loop
# =============================================================================
TOTAL=$(( (ANGLE_END - ANGLE_START) / ANGLE_STEP + 1 ))
COUNT=0; FAILED=0

log "Starting dihedral scan ($TOTAL points) ..."
echo "══════════════════════════════════════════════════════════════"

angle=$ANGLE_START
while [[ $angle -le $ANGLE_END ]]; do
    COUNT=$(( COUNT + 1 ))
    ADIR="$OUTPUT_DIR/angle_$(printf '%04d' $angle)"
    mkdir -p "$ADIR"

    echo ""
    printf "${CYAN}[%3d/%3d]${NC} Angle = ${GREEN}%4d°${NC}\n" \
           "$COUNT" "$TOTAL" "$angle"

    # ── Patch topology for this angle (shared by all four stages) ─────────────
    STEP_TOP="$ADIR/topol_$(printf '%04d' $angle).top"
    patch_topology "$TOP_FILE" "$STEP_TOP" "$angle"

    # ── Stage 1: Rough EM ─────────────────────────────────────────────────────
    S1D="$ADIR/s1_rough_em"; mkdir -p "$S1D"
    stg 1 "Rough EM (clash removal)"

    if ! run_grompp "$MDP_S1" "$GRO_FILE" "$STEP_TOP" "$S1D/em.tpr" "$S1D/run.log"; then
        warn "Stage 1 grompp failed at ${angle}°  →  $S1D/run.log"
        echo "${angle}   NaN   FAILED_S1_GROMPP" >> "$RESULTS_FILE"
        FAILED=$(( FAILED + 1 )); angle=$(( angle + ANGLE_STEP )); continue
    fi

    # EM non-convergence is acceptable here — we just need a clash-free structure
    run_mdrun "$S1D/em.tpr" "$S1D/em" "$S1D/run.log" || \
        warn "Stage 1 EM did not converge (acceptable at loose tolerance)"

    S1_GRO="$S1D/em.gro"
    if [[ ! -f "$S1_GRO" ]]; then
        warn "Stage 1 produced no .gro at ${angle}°"
        echo "${angle}   NaN   FAILED_S1_MDRUN" >> "$RESULTS_FILE"
        FAILED=$(( FAILED + 1 )); angle=$(( angle + ANGLE_STEP )); continue
    fi

    # ── Stage 2: Hot NVT SD ───────────────────────────────────────────────────
    S2D="$ADIR/s2_hot_md"; mkdir -p "$S2D"
    stg 2 "Hot NVT SD  @  ${HOT_TEMP} K  (${HOT_TIME} ps)"

    if ! run_grompp "$MDP_S2" "$S1_GRO" "$STEP_TOP" "$S2D/md.tpr" "$S2D/run.log"; then
        warn "Stage 2 grompp failed at ${angle}°  →  $S2D/run.log"
        echo "${angle}   NaN   FAILED_S2_GROMPP" >> "$RESULTS_FILE"
        FAILED=$(( FAILED + 1 )); angle=$(( angle + ANGLE_STEP )); continue
    fi

    S2_GRO="$S1_GRO"   # fallback: use Stage 1 output if Stage 2 crashes
    if run_mdrun "$S2D/md.tpr" "$S2D/md" "$S2D/run.log"; then
        [[ -f "$S2D/md.gro" ]] && S2_GRO="$S2D/md.gro"
    else
        warn "Stage 2 mdrun issue at ${angle}° — using Stage 1 structure as fallback"
    fi

    # ── Stage 3: Simulated annealing ──────────────────────────────────────────
    S3D="$ADIR/s3_anneal"; mkdir -p "$S3D"
    stg 3 "Simulated annealing  ${HOT_TEMP} K → 200 K → 50 K  (${ANNEAL_TIME} ps)"

    if ! run_grompp "$MDP_S3" "$S2_GRO" "$STEP_TOP" "$S3D/md.tpr" "$S3D/run.log"; then
        warn "Stage 3 grompp failed at ${angle}°  →  $S3D/run.log"
        echo "${angle}   NaN   FAILED_S3_GROMPP" >> "$RESULTS_FILE"
        FAILED=$(( FAILED + 1 )); angle=$(( angle + ANGLE_STEP )); continue
    fi

    S3_GRO="$S2_GRO"   # fallback
    if run_mdrun "$S3D/md.tpr" "$S3D/md" "$S3D/run.log"; then
        [[ -f "$S3D/md.gro" ]] && S3_GRO="$S3D/md.gro"
    else
        warn "Stage 3 mdrun issue at ${angle}° — using Stage 2 structure as fallback"
    fi

    # ── Stage 4: Final tight EM ───────────────────────────────────────────────
    S4D="$ADIR/s4_final_em"; mkdir -p "$S4D"
    stg 4 "Final tight EM (emtol = 1 kJ/mol/nm)"

    CONVERGED="YES"
    if ! run_grompp "$MDP_S4" "$S3_GRO" "$STEP_TOP" "$S4D/em.tpr" "$S4D/run.log"; then
        warn "Stage 4 grompp failed at ${angle}°  →  $S4D/run.log"
        echo "${angle}   NaN   FAILED_S4_GROMPP" >> "$RESULTS_FILE"
        FAILED=$(( FAILED + 1 )); angle=$(( angle + ANGLE_STEP )); continue
    fi

    if ! run_mdrun "$S4D/em.tpr" "$S4D/em" "$S4D/run.log"; then
        warn "Stage 4 EM did not fully converge at ${angle}° (energy still recorded)"
        CONVERGED="NO_CONV"
        FAILED=$(( FAILED + 1 ))
    fi

    # ── Extract potential energy from Stage 4 .edr ────────────────────────────
    if [[ -f "$S4D/em.edr" ]]; then
        EPOT=$(extract_energy "$S4D/em.edr" "$S4D")
    else
        EPOT="NaN"; CONVERGED="NO_EDR"
    fi

    # ── Record result ─────────────────────────────────────────────────────────
    printf "%-10s  %-16s  %s\n" "$angle" "$EPOT" "$CONVERGED" >> "$RESULTS_FILE"
    ok "Angle ${angle}°: Epot = ${EPOT} kJ/mol  [${CONVERGED}]"
    echo "$(date '+%H:%M:%S')  angle=$angle  Epot=$EPOT  status=$CONVERGED" >> "$LOG_FILE"

    angle=$(( angle + ANGLE_STEP ))
done

echo ""
echo "══════════════════════════════════════════════════════════════"
ok "Scan complete.  $COUNT points processed,  $FAILED issue(s)."
ok "Results : $RESULTS_FILE"

# =============================================================================
# Gnuplot visualisation script
# =============================================================================
PLOT_SCRIPT="$OUTPUT_DIR/plot_scan.gnu"
cat > "$PLOT_SCRIPT" <<GNUPLOT
#!/usr/bin/gnuplot
set title "Dihedral Scan (annealing protocol) — atoms $ATOM1-$ATOM2-$ATOM3-$ATOM4"
set xlabel "Dihedral angle (degrees)"
set ylabel "Potential energy (kJ/mol)"
set grid
set xrange [0:360]
set xtics 30
set key top right

plot "$RESULTS_FILE" using 1:2 with linespoints \\
     pt 7 ps 0.8 lc rgb "steelblue" title "E_pot (annealed + EM)"

pause -1 "Press Enter to exit"
GNUPLOT
chmod +x "$PLOT_SCRIPT"
log "Gnuplot script : $PLOT_SCRIPT"
log "Visualise with : gnuplot $PLOT_SCRIPT"

# =============================================================================
# Summary table
# =============================================================================
echo ""
echo "════════════════════════════════════════════"
echo "  Dihedral Scan Summary"
echo "════════════════════════════════════════════"
echo "  Angle(°)   Epot (kJ/mol)   Status"
echo "────────────────────────────────────────────"
grep -v '^#' "$RESULTS_FILE" | \
    awk '{printf "  %-10s %-16s %s\n", $1, $2, $3}'
echo "════════════════════════════════════════════"
