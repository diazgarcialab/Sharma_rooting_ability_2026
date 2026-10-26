#!/usr/bin/env bash
# =============================================================================
# gemma_rooting.sh -- GWAS for every adventitious-rooting trait
#
# RUN WITH bash, NOT sh:      bash gemma_rooting.sh
#
#   bash gemma_rooting.sh            resolve + prep + gwas, resuming as needed
#   bash gemma_rooting.sh resolve    pick one sequencing library per accession
#   bash gemma_rooting.sh prep       genotype file + kinship + PCA
#   bash gemma_rooting.sh gwas       the five association tests
#   bash gemma_rooting.sh clean      delete everything derived and start over
#
#   TRAITS="callus_stage root_weight" bash gemma_rooting.sh gwas
#   FORCE=1 bash gemma_rooting.sh    ignore existing outputs and redo them
# =============================================================================

set -euo pipefail

# ---- configuration ----------------------------------------------------------
WORKDIR=${WORKDIR:-/data6/diazgarcia/rooting}
VCF_NAME=${VCF_NAME:-VITVgdSC2_v1.0.Hap1.varonly.hardfil.miss0.25.maf0.01.bisnp.norep.g.vcf.gz}
# A copy in the working directory wins over one in the chloride folder.
if [[ -z "${VCF:-}" ]]; then
  if [[ -s "$WORKDIR/$VCF_NAME" ]]; then VCF="$WORKDIR/$VCF_NAME"
  else VCF="/data6/diazgarcia/chloride/$VCF_NAME"; fi
fi

THREADS=${THREADS:-28}
MAF=${MAF:-0.05}                 # the chloride run used 0.1
N_PCS=${N_PCS:-10}
PLINK_MEM=${PLINK_MEM:-}         # MB; empty means do not pass --memory at all
PLINK=${PLINK:-plink1.9}
GEMMA=${GEMMA:-gemma}
BCFTOOLS=${BCFTOOLS:-bcftools}
FORCE=${FORCE:-0}
VERIFY_BLOCKS=${VERIFY_BLOCKS:-0}   # 1 = run bgzip -t over the whole VCF first

# Passing `--memory ''` makes PLINK abort, so build the flag once and expand it
# with the ${arr[@]+...} guard that keeps `set -u` happy on an empty array.
PLINK_MEMFLAG=()
[[ -n "$PLINK_MEM" ]] && PLINK_MEMFLAG=(--memory "$PLINK_MEM")

STAGE=${1:-all}
cd "$WORKDIR"
LOGFILE="$WORKDIR/gemma_rooting.$(date +%Y%m%d_%H%M%S).log"

log()  { printf '\n[%s] %s\n' "$(date +%H:%M:%S)" "$*" | tee -a "$LOGFILE"; }
info() { printf '        %s\n' "$*" | tee -a "$LOGFILE"; }
die()  { printf '\nERROR: %s\n' "$*" | tee -a "$LOGFILE" >&2; exit 1; }
stale() { [[ "$FORCE" == "1" || ! -s "$1" ]]; }

DERIVED=(geno.bed geno.bim geno.fam geno.fam.orig geno.log geno.nosex
         geno_pca.eigenval geno_pca.eigenvec geno_pca.log geno_pca.nosex
         covariates_numeric.txt vcf_samples.txt accessions_not_in_vcf.txt
         individuals_blups.txt sample_ids.txt fam_phenotypes.txt
         filtered_phenotypes_blups.txt)

clean() {
  log "removing derived files"
  rm -f "${DERIVED[@]}" gemma_rooting.*.log 2>/dev/null || true
  rm -rf output 2>/dev/null || true
  info "kept: fam_phenotypes_all.txt, trait_index.tsv, accession_libraries.csv,"
  info "      pheno_*.txt and the VCF"
  info "run  bash gemma_rooting.sh  to start over"
}

if [[ "$STAGE" == "clean" ]]; then clean; exit 0; fi

# ---- preflight --------------------------------------------------------------
for exe in "$PLINK" "$GEMMA" awk; do
  command -v "$exe" >/dev/null 2>&1 || die "$exe not found on PATH"
done
for f in fam_phenotypes_all.txt trait_index.tsv; do
  [[ -s "$f" ]] || die "missing input: $WORKDIR/$f
       Copy 02_clean_data/gwas/ up from the laptop first."
done
[[ -s "$VCF" ]] || die "VCF not found: $VCF"

N_TRAITS=$(awk -F'\t' 'NR==1{print NF-4; exit}' fam_phenotypes_all.txt)  # minus id,id,acc,dummy
N_PHENO_COLS=$(( N_TRAITS + 1 ))                                   # traits + dummy
AVAIL_GB=$(df -Pk . | awk 'NR==2{printf "%d", $4/1024/1024}')

log "log file: $LOGFILE"
info "VCF          : $VCF"
info "               $(du -h "$VCF" | cut -f1), modified $(date -r "$VCF" '+%Y-%m-%d %H:%M')"
info "free space   : ${AVAIL_GB} GB"
info "libraries    : $(wc -l < fam_phenotypes_all.txt)"
info "accessions   : $(awk -F'\t' '{print $3}' fam_phenotypes_all.txt | sort -u | wc -l)"
info "traits       : $N_TRAITS  (+1 dummy column = $N_PHENO_COLS phenotype columns)"
[[ "$AVAIL_GB" -ge 10 ]] || die "less than 10 GB free; the .bed will not fit"

# ---- VCF integrity ----------------------------------------------------------
# A damaged container costs minutes to discover through PLINK and a second here.
check_vcf() {
  local magic eofblk
  magic=$(head -c 16 "$VCF" | od -An -tx1 | tr -d ' \n')
  [[ "$magic" == 1f8b* ]] || die "$VCF is not gzip-compressed"
  if [[ "$magic" == *4243* ]]; then
    eofblk=$(tail -c 28 "$VCF" | od -An -tx1 | tr -d ' \n')
    [[ "$eofblk" == "1f8b08040000000000ff0600424302001b0003000000000000000000" ]] \
      || die "$VCF is missing its BGZF EOF block -- the copy is incomplete.
       Copy it again from backup and re-run."
    info "BGZF container, EOF block present"
  fi
  if [[ "$VERIFY_BLOCKS" == "1" ]] && command -v bgzip >/dev/null 2>&1; then
    info "VERIFY_BLOCKS=1: bgzip -t over the whole file, a few minutes ..."
    bgzip -t "$VCF" 2>/tmp/bgzip_check.$$ \
      || die "$VCF has a damaged block: $(head -1 /tmp/bgzip_check.$$)
       Copy it again from backup."
    rm -f /tmp/bgzip_check.$$
    info "  all blocks decompress"
  fi
}

# ---- stage: resolve ---------------------------------------------------------
vcf_sample_list() {
  if command -v "$BCFTOOLS" >/dev/null 2>&1; then
    "$BCFTOOLS" query -l "$VCF" > vcf_samples.txt
  else
    # header-only fallback when bcftools is unavailable
    zcat "$VCF" 2>/dev/null | awk '/^#CHROM/ {
      for (i = 10; i <= NF; i++) print $i; exit }' > vcf_samples.txt
  fi
  [[ -s vcf_samples.txt ]] || die "could not read any sample names from $VCF"
}

resolve() {
  if ! stale individuals_blups.txt && ! stale fam_phenotypes.txt && \
     ! stale vcf_samples.txt; then
    info "already resolved ($(wc -l < individuals_blups.txt) accessions) -- skipping"
    info "  FORCE=1 to redo"
    return
  fi
  log "resolving one sequencing library per accession"
  vcf_sample_list
  info "samples in the VCF: $(wc -l < vcf_samples.txt)"

  # fam_phenotypes_all.txt is sorted by accession then library index, so the
  # first library present in the VCF is the lowest-numbered one.
  #
  # It is TAB-separated because 38 accession names contain a space, hence the
  # two-file awk below splits the second file on tabs (FS is set in the rule for
  # FNR == NR's file too, but the VCF sample list is one field per line so it
  # does not care). Output is space-delimited, which is what PLINK and GEMMA want.
  rm -f accessions_not_in_vcf.txt
  awk -F'\t' -v missing=accessions_not_in_vcf.txt '
    NR == FNR { inv[$1] = 1; next }
    {
      acc = $3
      seen[acc] = 1
      if (!(acc in taken) && ($1 in inv)) {
        taken[acc] = 1
        printf "%s %s", $1, $2
        for (i = 4; i <= NF; i++) printf " %s", $i
        printf "\n"
      }
    }
    END { for (a in seen) if (!(a in taken)) print a > missing }
  ' vcf_samples.txt fam_phenotypes_all.txt > fam_phenotypes.txt

  awk '{print $1, $2}' fam_phenotypes.txt > individuals_blups.txt
  awk '{print $2}'     fam_phenotypes.txt > sample_ids.txt
  # drop the dummy column for the plain phenotype file
  awk '{ printf "%s %s", $1, $2
         for (i = 4; i <= NF; i++) printf " %s", $i
         printf "\n" }' fam_phenotypes.txt > filtered_phenotypes_blups.txt

  local n_res n_acc n_missing
  n_res=$(wc -l < individuals_blups.txt)
  n_acc=$(awk -F'\t' '{print $3}' fam_phenotypes_all.txt | sort -u | wc -l)
  n_missing=0
  [[ -s accessions_not_in_vcf.txt ]] && n_missing=$(wc -l < accessions_not_in_vcf.txt)

  info "accessions with a library in the VCF : $n_res / $n_acc"
  if [[ "$n_missing" -gt 0 ]]; then
    sort -o accessions_not_in_vcf.txt accessions_not_in_vcf.txt
    info "accessions with NO library in the VCF: $n_missing"
    info "  listed in $WORKDIR/accessions_not_in_vcf.txt"
    head -20 accessions_not_in_vcf.txt | sed 's/^/          /' | tee -a "$LOGFILE"
    [[ "$n_missing" -gt 20 ]] && info "          ... and $((n_missing - 20)) more"
    info "  The GWAS runs on $n_res accessions. Report that number, not $n_acc."
  else
    rm -f accessions_not_in_vcf.txt
    info "every accession has a library in the VCF"
  fi
}

# ---- stage: prep ------------------------------------------------------------
build_geno() {
  if ! stale geno.bed && ! stale geno.bim && [[ -s geno.fam.orig ]]; then
    info "geno.bed present ($(wc -l < geno.bim) variants) -- skipping the import"
    return
  fi
  check_vcf
  log "importing the VCF with PLINK (the slow step)"
  "$PLINK" --vcf "$VCF" \
           --keep individuals_blups.txt \
           --make-bed --out geno \
           --threads "$THREADS" ${PLINK_MEMFLAG[@]+"${PLINK_MEMFLAG[@]}"} \
           --allow-extra-chr --maf "$MAF"
  [[ -s geno.bed ]] || die "PLINK did not produce geno.bed"
  cp geno.fam geno.fam.orig
  local n_want n_got
  n_want=$(wc -l < individuals_blups.txt); n_got=$(wc -l < geno.fam)
  info "geno: $(wc -l < geno.bim) variants, $n_got samples"
  [[ "$n_got" -eq "$n_want" ]] || info "  NOTE: asked for $n_want, kept $n_got"
}

patch_fam() {
  [[ -s geno.fam.orig ]] || die "geno.fam.orig missing -- rerun with FORCE=1"
  log "writing the phenotype columns into geno.fam"
  # Match on IID and keep the .fam row ORDER: GEMMA pairs the .fam, the kinship
  # matrix and the covariate file by row, not by ID.
  awk -v nt="$N_PHENO_COLS" '
    NR == FNR { s = ""; for (i = 3; i <= NF; i++) s = s " " $i; ph[$2] = s; next }
    {
      printf "%s %s %s %s %s", $1, $2, $3, $4, $5
      if ($2 in ph) printf "%s", ph[$2]
      else for (i = 1; i <= nt; i++) printf " NA"
      printf "\n"
    }' fam_phenotypes.txt geno.fam.orig > geno.fam.tmp
  mv geno.fam.tmp geno.fam
  local n_ph; n_ph=$(awk '$6 != "NA"' geno.fam | wc -l)
  info "geno.fam: $(wc -l < geno.fam) samples, $n_ph carry phenotypes"
  [[ "$n_ph" -gt 0 ]] || die "no sample in geno.fam matched fam_phenotypes.txt"
}

build_kinship() {
  if ! stale output/kinship.sXX.txt; then
    info "output/kinship.sXX.txt present -- skipping"; return
  fi
  log "kinship matrix (standardised, -gk 2)"
  "$GEMMA" -bfile geno -gk 2 -n 1 -o kinship     # -n 1 is the complete dummy
  [[ -s output/kinship.sXX.txt ]] || die "GEMMA did not produce the kinship matrix"
}

build_pca() {
  if stale geno_pca.eigenvec; then
    log "principal components"
    "$PLINK" --bfile geno --pca "$N_PCS" --out geno_pca \
             --allow-extra-chr ${PLINK_MEMFLAG[@]+"${PLINK_MEMFLAG[@]}"}
  else
    info "geno_pca.eigenvec present -- skipping"
  fi
  log "covariate file (intercept + $N_PCS PCs, in .fam order)"
  awk 'NR == FNR { s = ""; for (i = 3; i <= NF; i++) s = s " " $i; pc[$2] = s; next }
       { if (!($2 in pc)) { print "missing PCs for " $2 > "/dev/stderr"; exit 1 }
         print "1" pc[$2] }' geno_pca.eigenvec geno.fam > covariates_numeric.txt
  [[ "$(wc -l < covariates_numeric.txt)" -eq "$(wc -l < geno.fam)" ]] \
    || die "covariate rows do not match geno.fam"
  info "covariates: $(awk 'NR==1{print NF}' covariates_numeric.txt) columns x $(wc -l < geno.fam) rows"
}

prep() { build_geno; patch_fam; build_kinship; build_pca; log "prep done"; }

# ---- stage: gwas ------------------------------------------------------------
gwas() {
  [[ -s output/kinship.sXX.txt ]] || die "run the prep stage first"
  [[ -s covariates_numeric.txt ]] || die "covariates_numeric.txt missing; run prep"
  local wanted=${TRAITS:-}

  while IFS=$'\t' read -r n_index trait slug pheno_file n_samples value_type; do
    [[ "$n_index" == "n_index" ]] && continue
    if [[ -n "$wanted" && " $wanted " != *" $slug "* ]]; then continue; fi

    local out="gwas_${slug}_blup_PCA${N_PCS}"
    if ! stale "output/${out}.assoc.txt"; then
      info "$trait done ($(( $(wc -l < "output/${out}.assoc.txt") - 1 )) SNPs) -- skipping"
      continue
    fi
    log "GWAS: $trait   (-n $n_index)"
    "$GEMMA" -bfile geno -k output/kinship.sXX.txt -lmm 1 \
             -c covariates_numeric.txt -n "$n_index" -o "$out"
    [[ -s "output/${out}.log.txt" ]] && grep -E \
      'number of analyzed individuals|number of analyzed SNPs|pve estimate|se\(pve\)' \
      "output/${out}.log.txt" | sed 's/^/        /' | tee -a "$LOGFILE" || true
  done < trait_index.tsv
  log "association tests done"
}

case "$STAGE" in
  resolve) resolve ;;
  prep)    prep ;;
  gwas)    gwas ;;
  all)     resolve; prep; gwas ;;
  *)       die "unknown stage '$STAGE' (use resolve, prep, gwas, clean or all)" ;;
esac

log "results in $WORKDIR/output/"
ls -1 "$WORKDIR"/output/gwas_*.assoc.txt 2>/dev/null | sed 's/^/        /' || true
