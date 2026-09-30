#!/usr/bin/env bash
# Quick sanity check that every tool in the lab is on the PATH.
tools=(nextflow nf-core nf-test snakemake conda mamba docker java
       fastqc fastp multiqc seqkit bwa minimap2 salmon
       samtools bcftools bgzip tabix bedtools python)

missing=0
echo "Checking tools:"
for t in "${tools[@]}"; do
  if command -v "$t" > /dev/null 2>&1; then
    printf "  OK       %s\n" "$t"
  else
    printf "  MISSING  %s\n" "$t"
    missing=$((missing + 1))
  fi
done

if [ "$missing" -eq 0 ]; then
  echo "All tools found."
else
  echo "$missing tool(s) missing. Try: conda activate bioinfo"
  exit 1
fi
