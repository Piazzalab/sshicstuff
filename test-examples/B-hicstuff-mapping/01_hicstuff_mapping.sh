#!/usr/bin/env bash
set -Eeuo pipefail

# ------------------------------------------------------------------------------
# Hi-C processing pipeline
# Uses: parasplit, hicstuff, cooler, bowtie2-build
#
# Steps:
#   1. Optional read digestion with parasplit
#   2. hicstuff pipeline -> graal matrix
#   3. rebin at 1 kb
#   4. convert to cooler
#   5. ICE balance
#   6. zoomify to mcool
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# User parameters
# ------------------------------------------------------------------------------

thread=4
mode="normal"                  # parasplit | cutsite | normal
enzymes="DpnII,HinfI"
quality=20

fastqdir="./fastq"
genome_fasta="./S288c_DSB_LY_Capture_artificial.fa"   # <-- your genome FASTA
bt2_index_dir="./bt2_index"            # where Bowtie2 index will be stored
outputdir="./localtmp"
outdir_digest="${fastqdir}"

SAMP=("AD433_sub4M")

r1ext=".end1"
r2ext=".end2"
fqext=".fastq.gz"

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

ts() {
    date '+%Y-%m-%d %H:%M:%S'
}

log() {
    echo "[$(ts)] $*"
}

die() {
    echo "[$(ts)] ERROR: $*" >&2
    exit 1
}

check_dependency() {
    command -v "$1" >/dev/null 2>&1 || die "$1 not found in PATH."
}

# ------------------------------------------------------------------------------
# Dependency checks
# ------------------------------------------------------------------------------

check_dependency hicstuff
check_dependency cooler
check_dependency parasplit
check_dependency bowtie2-build

# ------------------------------------------------------------------------------
# Input checks
# ------------------------------------------------------------------------------

[[ -f "${genome_fasta}" ]] || die "Genome FASTA not found: ${genome_fasta}"
[[ "${mode}" == "parasplit" || "${mode}" == "cutsite" || "${mode}" == "normal" ]] \
    || die "Invalid mode: ${mode}"

for sample in "${SAMP[@]}"; do
    in_r1="${fastqdir}/${sample}${r1ext}${fqext}"
    in_r2="${fastqdir}/${sample}${r2ext}${fqext}"

    [[ -f "${in_r1}" ]] || die "R1 file not found: ${in_r1}"
    [[ -f "${in_r2}" ]] || die "R2 file not found: ${in_r2}"
done

mkdir -p "${outputdir}" "${outdir_digest}" "${bt2_index_dir}"

# ------------------------------------------------------------------------------
# Genome / Bowtie2 index setup
# ------------------------------------------------------------------------------

genome_base="$(basename "${genome_fasta}")"      # S288c_DSB_edited.fa
genome_name="${genome_base%.*}"                  # S288c_DSB_edited
bt2_prefix="${bt2_index_dir}/${genome_name}"

# Build Bowtie2 index only if missing
if [[ ! -e "${bt2_prefix}.1.bt2" && ! -e "${bt2_prefix}.1.bt2l" ]]; then
    log "Bowtie2 index not found. Building index from ${genome_fasta}"
    bowtie2-build "${genome_fasta}" "${bt2_prefix}"
else
    log "Bowtie2 index found: ${bt2_prefix}"
fi

# hicstuff version string for filenames
hcsver_raw="$(hicstuff -v 2>/dev/null || hicstuff --version 2>/dev/null || echo unknown)"
hcsver="v${hcsver_raw}"
hcsver="${hcsver// /_}"
hcsver="${hcsver//./}"
hcsver="${hcsver//\//_}"

# ------------------------------------------------------------------------------
# Run
# ------------------------------------------------------------------------------

log "Started reads_to_hic.sh"

echo "Parameters:"
for v in thread mode enzymes quality fastqdir genome_fasta bt2_index_dir outputdir outdir_digest; do
    declare -n ref="$v"
    echo "  [$v] ${ref}"
done

echo "Tool versions:"
for t in parasplit hicstuff cooler bowtie2-build; do
    ver="$($t --version 2>/dev/null | head -n 1 || true)"
    echo "  [$t] ${ver}"
done

echo "Samples:"
printf '  %s\n' "${SAMP[@]}"

for sample in "${SAMP[@]}"; do
    log "Started sample: ${sample}"

    sdir="${outputdir}/${sample}"
    scdir="${sdir}/Cool"
    logdir="${sdir}/logs"
    mkdir -p "${sdir}" "${scdir}" "${logdir}"

    sample_mode="${mode}"   # local copy, do not mutate global mode

    out_px="${sample}_${genome_name}_${sample_mode}_q${quality}_${hcsver}"

    in_r1="${fastqdir}/${sample}${r1ext}${fqext}"
    in_r2="${fastqdir}/${sample}${r2ext}${fqext}"

    # --------------------------------------------------------------------------
    # Optional parasplit digestion
    # --------------------------------------------------------------------------
    if [[ "${sample_mode}" == "parasplit" ]]; then
        log "  Digesting reads with parasplit"

        para_r1="${outdir_digest}/${sample}_parasplit${r1ext}${fqext}"
        para_r2="${outdir_digest}/${sample}_parasplit${r2ext}${fqext}"

        parasplit \
            --source_forward="${in_r1}" \
            --source_reverse="${in_r2}" \
            --output_forward="${para_r1}" \
            --output_reverse="${para_r2}" \
            --list_enzyme="${enzymes}" \
            --mode="all" \
            --num_threads="${thread}" \
            > "${logdir}/parasplit.stdout.log" \
            2> "${logdir}/parasplit.stderr.log"

        in_r1="${para_r1}"
        in_r2="${para_r2}"

        # After digestion, hicstuff mapping mode should be normal
        sample_mode="normal"
    fi

    # --------------------------------------------------------------------------
    # hicstuff pipeline
    # --------------------------------------------------------------------------
    log "  Running hicstuff pipeline"

    hicstuff pipeline \
        --duplicates \
        --distance-law \
        --filter \
        --force \
        --no-cleanup \
        --plot \
        --matfmt=graal \
        --enzyme="${enzymes}" \
        --mapping="${sample_mode}" \
        --quality-min="${quality}" \
        --threads="${thread}" \
        --outdir "${sdir}" \
        --genome="${bt2_prefix}" \
        "${in_r1}" "${in_r2}" \
        > "${logdir}/hicstuff_pipeline.stdout.log" \
        2> "${logdir}/hicstuff_pipeline.stderr.log"

    # Expected hicstuff outputs (native hicstuff names)
    raw_frag_contacts="${sdir}/abs_fragments_contacts_weighted.txt"
    raw_frag_list="${sdir}/fragments_list.txt"
    raw_contigs="${sdir}/info_contigs.txt"

    [[ -f "${raw_frag_contacts}" ]] || die "Missing file after hicstuff pipeline: ${raw_frag_contacts}"
    [[ -f "${raw_frag_list}" ]] || die "Missing file after hicstuff pipeline: ${raw_frag_list}"
    [[ -f "${raw_contigs}" ]] || die "Missing file after hicstuff pipeline: ${raw_contigs}"

    # Rename outputs to sample-specific filenames for clarity
    frag_contacts="${sdir}/${sample}_abs_graal_fragments_weighted.txt"
    frag_list="${sdir}/${sample}_digested_fragments_list.txt"
    contigs="${sdir}/${sample}_info_contigs.txt"

    mv -f "${raw_frag_contacts}" "${frag_contacts}"
    mv -f "${raw_frag_list}" "${frag_list}"
    mv -f "${raw_contigs}" "${contigs}"
    # --------------------------------------------------------------------------
    # Rebin at 1 kb
    # --------------------------------------------------------------------------
    log "  Rebinning matrix at 1 kb"

    rebinned_prefix="${sdir}/${out_px}_1kb"

    hicstuff rebin \
        --binning=1kb \
        --force \
        --frags="${frag_list}" \
        --chroms="${contigs}" \
        "${frag_contacts}" \
        "${rebinned_prefix}" \
        > "${logdir}/hicstuff_rebin.stdout.log" \
        2> "${logdir}/hicstuff_rebin.stderr.log"

    [[ -f "${rebinned_prefix}.mat.tsv" ]] || die "Missing rebinned matrix: ${rebinned_prefix}.mat.tsv"
    [[ -f "${rebinned_prefix}.frags.tsv" ]] || die "Missing rebinned fragments: ${rebinned_prefix}.frags.tsv"
    [[ -f "${rebinned_prefix}.chr.tsv" ]] || die "Missing rebinned chroms: ${rebinned_prefix}.chr.tsv"

    # --------------------------------------------------------------------------
    # Convert to cooler
    # --------------------------------------------------------------------------
    log "  Converting rebinned matrix to cooler"

    hicstuff convert \
        --force \
        --frags="${rebinned_prefix}.frags.tsv" \
        --chroms="${rebinned_prefix}.chr.tsv" \
        "${rebinned_prefix}.mat.tsv" \
        "${scdir}/${out_px}_1kb" \
        > "${logdir}/hicstuff_convert.stdout.log" \
        2> "${logdir}/hicstuff_convert.stderr.log"

    cool_1kb="${scdir}/${out_px}_1kb.cool"
    [[ -f "${cool_1kb}" ]] || die "Missing cooler file: ${cool_1kb}"

    # --------------------------------------------------------------------------
    # Balance
    # --------------------------------------------------------------------------
    log "  ICE-balancing matrix"

    cooler balance \
        --nproc "${thread}" \
        "${cool_1kb}" \
        > "${logdir}/cooler_balance.stdout.log" \
        2> "${logdir}/cooler_balance.stderr.log"

    # --------------------------------------------------------------------------
    # Zoomify
    # --------------------------------------------------------------------------
    log "  Zoomifying matrix"

    mcool="${scdir}/${out_px}.mcool"

    cooler zoomify \
        --balance \
        --nproc "${thread}" \
        --resolutions 1000,2000,5000,10000,15000,20000 \
        --out "${mcool}" \
        "${cool_1kb}" \
        > "${logdir}/cooler_zoomify.stdout.log" \
        2> "${logdir}/cooler_zoomify.stderr.log"

    [[ -f "${mcool}" ]] || die "Missing mcool file: ${mcool}"

    log "Finished sample: ${sample}"
done

log "Done."
