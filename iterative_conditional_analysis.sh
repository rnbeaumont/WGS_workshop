#!/bin/bash

set -euo pipefail

usage() {
  cat <<'EOF'
Usage: step3_iterative_conditional_analysis.sh <sig_thresh> <initial_regenie_file> <seed_condition_list> <mask_file> <set_list_file> <annotation_file> [options]

Required positional arguments:
  sig_thresh             LOG10P threshold used to select significant aggregates
  initial_regenie_file    Starting .regenie output file for the iterative loop
  seed_condition_list     File of single-variant hits to seed conditioning
  mask_file               REGENIE mask definition file
  set_list_file           REGENIE set-list file
  annotation_file        REGENIE annotation file

Options:
  --bed-prefix PATH       PLINK bed prefix to use for the extract step and iterative reruns
  --afreq-file PATH       PLINK .afreq file used to derive rare-variant MAF filtering
  --pheno-file PATH       Phenotype file used by REGENIE
  --pheno-col NAME        Phenotype column name used by REGENIE
  --covar-file PATH       Covariate file used by REGENIE
  --cat-covar-list LIST   Comma-separated categorical covariates for REGENIE
  --pred-file PATH        REGENIE step 1 prediction list
  --out-prefix PATH       Prefix for outputs from the initial extracted PLINK bed
  --threads N             Threads passed to REGENIE
  --bsize N               Block size passed to REGENIE
  --max-condition-vars N   Maximum condition variables passed to REGENIE
  --vc-tests VALUE        REGENIE vc-tests value
  --vc-maxAAF VALUE       REGENIE vc-maxAAF value
  --joint VALUE           REGENIE joint test value
  --aaf-bins VALUE        REGENIE aaf-bins value
  --chr VALUE             Chromosome passed to REGENIE
  --conda-env NAME        Conda environment used with conda run fallback (default: regenie)
  --regenie-bin CMD       REGENIE executable name or absolute path (default: regenie)
  --output_folder FOLDER  Output folder for PLINK files and REGENIE outputs from iterative runs
EOF
}

if [[ $# -lt 6 ]]; then
  usage >&2
  exit 1
fi

sig_thresh="$1"
initial_regenie_file="$2"
seed_condition_list="$3"
mask_file="$4"
set_list_file="$5"
annotation_file="$6"
shift 6

bed_prefix="ALL.chr2.build38"
pheno_file="simulated_phenotype.phenotype.tsv"
pheno_col="phenotype_sim"
covar_file="simulated_phenotype.covariate.tsv"
cat_covar_list="sex,ancestry"
out_prefix="${bed_prefix}_sig_aggregates_plus_condition"
threads=10
bsize=400
max_condition_vars=20000
vc_tests="acato-full"
vc_maxAAF=0.01
joint_test="acat"
aaf_bins=0.01
chr_value=2
conda_env_name="regenie"
regenie_bin="regenie"
output_folder="."
afreq_file="$output_folder/1kg_chr2.afreq"
pred_file="$output_folder/regenie_step1_pred.list"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --bed-prefix)
      bed_prefix="$2"
      shift 2
      ;;
    --afreq-file)
      afreq_file="$2"
      shift 2
      ;;
    --pheno-file)
      pheno_file="$2"
      shift 2
      ;;
    --pheno-col)
      pheno_col="$2"
      shift 2
      ;;
    --covar-file)
      covar_file="$2"
      shift 2
      ;;
    --cat-covar-list)
      cat_covar_list="$2"
      shift 2
      ;;
    --pred-file)
      pred_file="$2"
      shift 2
      ;;
    --out-prefix)
      out_prefix="$2"
      shift 2
      ;;
    --threads)
      threads="$2"
      shift 2
      ;;
    --bsize)
      bsize="$2"
      shift 2
      ;;
    --max-condition-vars)
      max_condition_vars="$2"
      shift 2
      ;;
    --vc-tests)
      vc_tests="$2"
      shift 2
      ;;
    --vc-maxAAF)
      vc_maxAAF="$2"
      shift 2
      ;;
    --joint)
      joint_test="$2"
      shift 2
      ;;
    --aaf-bins)
      aaf_bins="$2"
      shift 2
      ;;
    --chr)
      chr_value="$2"
      shift 2
      ;;
    --conda-env)
      conda_env_name="$2"
      shift 2
      ;;
    --regenie-bin)
      regenie_bin="$2"
      shift 2
      ;;
    --output_folder)
      output_folder="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

# Resolve REGENIE launcher for non-interactive shells:
# use executable on PATH when available, otherwise run inside conda env.
declare -a regenie_cmd
if command -v "$regenie_bin" >/dev/null 2>&1; then
  regenie_cmd=("$regenie_bin")
elif command -v conda >/dev/null 2>&1; then
  regenie_cmd=("conda" "run" "--no-capture-output" "-n" "$conda_env_name" "$regenie_bin")
else
  echo "Error: could not find '$regenie_bin' on PATH and 'conda' is unavailable for fallback." >&2
  exit 1
fi

echo "Using REGENIE command: ${regenie_cmd[*]}"

mkdir -p $output_folder

current_regenie_output="$initial_regenie_file"
log_file="$output_folder/conditional_analysis_iteration_log.txt"
selected_aggregates_file="$output_folder/final_selected_aggregates_by_iteration.tsv"

awk -v thresh="$sig_thresh" '
  NR>2{
    p=$12
    lp=tolower(p)
    if (p != "" && lp != "na" && lp != "nan" && lp != "-nan" && p+0==p && p>thresh) {
      print $3
    }
  }
' "$current_regenie_output" | tr '.' '\t' | cut -f 1-2 | sort | uniq | awk 'NF>1' > $output_folder/aggregates_significant_hits

awk '
  BEGIN{FS="[ \t]+"}
  NR==FNR{
    if(NF>=2){keep[$1"."$2]=1}
    next
  }
  ($2"."$3 in keep){
    print $1
  }
' $output_folder/aggregates_significant_hits "$annotation_file" > $output_folder/variants_from_significant_aggregates

cat $output_folder/variants_from_significant_aggregates "$seed_condition_list" | awk 'NF>0' | sort -u > $output_folder/variants_for_conditional_plink_extract

./plink2 \
  --bfile "$bed_prefix" \
  --extract $output_folder/variants_for_conditional_plink_extract \
  --make-bed \
  --out "$out_prefix"

awk '
  NR>1 && $5!="NA"{
    af=$5+0
    maf=(af<=0.5 ? af : 1-af)
    if(maf<0.011){print $2}
  }
' "$afreq_file" | sort -u > $output_folder/variants_maf_lt_0.011

cp "$seed_condition_list" $output_folder/condition_list_iter
last_completed_regenie_output="$current_regenie_output"
last_nonempty_sig_output="$current_regenie_output"
echo -e "iter\tsig_aggregates\taggregate_variants\trare_variants_added\tcondition_list_size\tregenie_output" > "$log_file"
echo -e "iter\tselected_aggregate\tselected_LOG10P\tsig_aggregates_in_source\tsource_regenie_output" > "$selected_aggregates_file"
echo "Starting iterative conditional analysis (sig_thresh=$sig_thresh)"
iter=1

while true; do
  current_sig_hits_all="$output_folder/aggregates_significant_hits_all_iter_${iter}"
  current_top_aggregate_info="$output_folder/top_aggregate_info_iter_${iter}.tsv"
  current_sig_hits="$output_folder/aggregates_significant_hits_iter_${iter}"
  current_sig_variants="$output_folder/variants_from_significant_aggregates_iter_${iter}"
  current_sig_rare_variants="$output_folder/rare_variants_from_significant_aggregates_iter_${iter}"

  awk -v thresh="$sig_thresh" '
    NR>2{
      p=$12
      lp=tolower(p)
      if (p != "" && lp != "na" && lp != "nan" && lp != "-nan" && p+0==p && p>thresh) {
        print $3
      }
    }
  ' "$current_regenie_output" | sort -u > "$current_sig_hits_all"

  if [[ ! -s "$current_sig_hits_all" ]]; then
    echo "Iteration ${iter}: no significant aggregates in $current_regenie_output; stopping."
    break
  fi

  n_sig_aggregates=$(wc -l < "$current_sig_hits_all")
  last_nonempty_sig_output="$current_regenie_output"

  awk -v thresh="$sig_thresh" '
    NR>2{
      p=$12
      lp=tolower(p)
      if (p != "" && lp != "na" && lp != "nan" && lp != "-nan" && p+0==p && p>thresh) {
        print $3"\t"$12
      }
    }
  ' "$current_regenie_output" | sort -k2,2gr | sed -n '1p' > "$current_top_aggregate_info"

  awk '{split($1,a,"."); if(length(a[1])>0 && length(a[2])>0){print a[1]"\t"a[2]}}' "$current_top_aggregate_info" > "$current_sig_hits"

  if [[ ! -s "$current_sig_hits" ]]; then
    echo "Iteration ${iter}: unable to parse top aggregate ID from $current_regenie_output; stopping."
    break
  fi

  selected_aggregate=$(awk 'NR==1{print $1}' "$current_top_aggregate_info")
  selected_log10p=$(awk 'NR==1{print $2}' "$current_top_aggregate_info")
  echo -e "${iter}\t${selected_aggregate}\t${selected_log10p}\t${n_sig_aggregates}\t${current_regenie_output}" >> "$selected_aggregates_file"

  awk '
    BEGIN{FS="[ \t]+"}
    NR==FNR{
      if(NF>=2){keep[$1"."$2]=1}
      next
    }
    ($2"."$3 in keep){
      print $1
    }
  ' "$current_sig_hits" "$annotation_file" | sort -u > "$current_sig_variants"

  awk '
    NR==FNR{rare[$1]=1; next}
    ($1 in rare){print $1}
  ' $output_folder/variants_maf_lt_0.011 "$current_sig_variants" > "$current_sig_rare_variants"

  if [[ ! -s "$current_sig_rare_variants" ]]; then
    n_agg_variants=$(wc -l < "$current_sig_variants")
    echo -e "${iter}\t${n_sig_aggregates}\t${n_agg_variants}\t0\t$(wc -l < $output_folder/condition_list_iter)\t$current_regenie_output" >> "$log_file"
    echo "Iteration ${iter}: significant aggregates found, but no variants with MAF < 0.011; stopping."
    break
  fi

  n_agg_variants=$(wc -l < "$current_sig_variants")
  n_rare_variants=$(wc -l < "$current_sig_rare_variants")

  cat $output_folder/condition_list_iter "$current_sig_rare_variants" | sort -u > $output_folder/condition_list_iter_next

  if cmp -s $output_folder/condition_list_iter $output_folder/condition_list_iter_next; then
    echo -e "${iter}\t${n_sig_aggregates}\t${n_agg_variants}\t0\t$(wc -l < $output_folder/condition_list_iter)\t$current_regenie_output" >> "$log_file"
    echo "Iteration ${iter}: no new conditional variants added; stopping."
    rm -f $output_folder/condition_list_iter_next
    break
  fi

  mv $output_folder/condition_list_iter_next $output_folder/condition_list_iter
  n_condition_vars=$(wc -l < $output_folder/condition_list_iter)

  echo -e "${iter}\t${n_sig_aggregates}\t${n_agg_variants}\t${n_rare_variants}\t${n_condition_vars}\t$current_regenie_output" >> "$log_file"
  echo "Iteration ${iter}: ${n_sig_aggregates} significant aggregates, ${n_rare_variants} rare variants added, condition list size now ${n_condition_vars}."

  n_bed_variants=$(wc -l < "${out_prefix}.bim")
  n_condition_in_bed=$(awk '
    NR==FNR{cond[$1]=1; next}
    ($2 in cond){n++}
    END{print n+0}
  ' $output_folder/condition_list_iter "${out_prefix}.bim")
  if [[ "$n_condition_in_bed" -ge "$n_bed_variants" ]]; then
    echo -e "${iter}\t${n_sig_aggregates}\t${n_agg_variants}\t0\t${n_condition_vars}\tSTOP_ALL_BED_VARIANTS_CONDITIONED" >> "$log_file"
    echo "Iteration ${iter}: condition list covers all ${n_bed_variants} variants in ${out_prefix}; stopping before REGENIE."
    break
  fi

  iter_out_prefix="$output_folder/regenie_step2_aggregates_cond_iter${iter}"
  if "${regenie_cmd[@]}" \
    --step 2 \
    --bed "$out_prefix" \
    --phenoFile "$pheno_file" \
    --phenoCol "$pheno_col" \
    --covarFile "$covar_file" \
    --catCovarList "$cat_covar_list" \
    --out "$iter_out_prefix" \
    --chr "$chr_value" \
    --pred "$pred_file" \
    --bsize "$bsize" \
    --apply-rint \
    --anno-file "$annotation_file" \
    --set-list "$set_list_file" \
    --mask-def "$mask_file" \
    --vc-tests "$vc_tests" \
    --vc-maxAAF "$vc_maxAAF" \
    --joint "$joint_test" \
    --aaf-bins "$aaf_bins" \
    --threads "$threads" \
    --check-burden-files \
    --max-condition-vars "$max_condition_vars" \
    --condition-list $output_folder/condition_list_iter; then
    :
  else
    regenie_exit_code=$?
    echo -e "${iter}\t${n_sig_aggregates}\t${n_agg_variants}\t0\t${n_condition_vars}\tREGENIE_FAILED_EXIT_${regenie_exit_code}" >> "$log_file"
    echo "Iteration ${iter}: REGENIE failed (exit code ${regenie_exit_code}); stopping iterative conditional analysis."
    break
  fi

  current_regenie_output="${iter_out_prefix}_phenotype_sim.regenie"
  if [[ ! -s "$current_regenie_output" ]]; then
    echo -e "${iter}\t${n_sig_aggregates}\t${n_agg_variants}\t0\t${n_condition_vars}\tREGENIE_OUTPUT_MISSING" >> "$log_file"
    echo "Iteration ${iter}: expected output $current_regenie_output was not created; stopping iterative conditional analysis."
    break
  fi

  last_completed_regenie_output="$current_regenie_output"
  iter=$((iter+1))
done

echo "Final selected aggregate list written to $selected_aggregates_file"
