#!/bin/bash
# Local only: no farm jobs. Existing products are never overwritten.
set -euo pipefail
[[ $# == 2 ]] || { echo 'Usage: bash tests/run_femto_lambda_validation.sh 13p5|3p85 OUTPUT_DIRECTORY' >&2; exit 2; }
study_energy="$1"
case "$study_energy" in 13p5|3p85) ;; *) exit 2;; esac
study_project=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)
cd "$study_project"
study_directory="$2"
[[ "$study_directory" == rootfile/femto_lambda_kf_validation_20260930 ]] || { echo 'Choose the dedicated validation destination in this frozen test script.' >&2; exit 2; }
study_input="config/picoDstList/auau${study_energy}_femtoLambda_validation_20260930.list"
for study_species in d t 3He 4He legacy standalone; do
  [[ ! -e "$study_directory/$study_energy/$study_species.root" ]] || { echo "Existing output $study_species" >&2; exit 2; }
done
for study_species in d t 3He 4He; do
  study_conf="config/mainconf/main_auau${study_energy}_anaFemtoLambda_${study_species}_KFParticle_highpurity.yaml"
  bash "script/singularity_run_anaFemtoLambda_${study_species}.sh" "$study_conf" "$study_input" \
    "$study_directory/$study_energy/$study_species.root" "validation_${study_energy}_${study_species}" 10000 \
    > "$study_directory/provenance/${study_energy}_${study_species}.log" 2>&1
  echo "COMPLETED new $study_energy $study_species"
done
study_legacy_conf="config/mainconf/main_auau${study_energy}_anaLambdaNuclearId.yaml"
bash tests/run_femto_lambda_root5.sh "$study_legacy_conf" \
  "analysis/run_anaLambdaNuclearId.C(\"$study_input\",\"$study_directory/$study_energy/legacy.root\",\"legacy_${study_energy}\",10000,\"$study_legacy_conf\")" \
  > "$study_directory/provenance/${study_energy}_legacy.log" 2>&1
rg -q 'StLambdaMaker::Finish\(\) processed 10000 events' "$study_directory/provenance/${study_energy}_legacy.log"
echo "COMPLETED original anaLambdaNuclearId $study_energy"
study_kf_conf="config/mainconf/main_auau${study_energy}_KFParticle_femto_control_20260930.yaml"
bash tests/run_femto_lambda_root5.sh "$study_kf_conf" \
  "analysis/run_anaLambda_KFParticle.C(\"$study_input\",\"$study_directory/$study_energy/standalone.root\",\"control_${study_energy}\",10000,\"$study_kf_conf\")" \
  > "$study_directory/provenance/${study_energy}_standalone.log" 2>&1
echo "COMPLETED standalone KF control $study_energy"
