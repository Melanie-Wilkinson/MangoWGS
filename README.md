# MangoWGS
225 Mangifera indica, raw vcf file with non-variants

1. SNP filtering from raw vcf
2. Population genetics analyses
3. Identifying deleterious mutations using SIFT

## 1. SNP filtering
### QF1
minGQ, Minimum genotype quality of 20 

```
vcftools --gzvcf ./qf4/1_rawVCF/4f90968ba585820a0142bd5dc0d5562d_raw_variants_1.vcf.gz --minGQ 20 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/1_GQ.recode.vcf.gz
```

minDP, Minimum depth per sample of 5
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/1_GQ.recode.vcf.gz --minDP 5 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/2_mindp.vcf.gz
```

Missing data per site: initially a relaxed threshold (i.e. 50%)
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/2_mindp.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/3_max-missing.vcf.gz
```

Create a list of missing data per individual:
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/3_max-missing.vcf.gz --missing-indv --out ./qf4/2_quality_filtering/4_
```
No individuals were removed as all have more than 85% of data
  
Maximum mean depth for removal of paralogues
1. Create a list of mean depth per site:
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/3_max-missing.vcf.gz --site-mean-depth --out ./qf4/2_quality_filtering/5_mean_depth
```
2. Filter for maximum mean depth of 50 for all sites
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/3_max-missing.vcf.gz --max-meanDP 50 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/5b_max_mdp.vcf.gz
```	


More stringent missing data per locus, filter for overall 20% missing data:
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/5b_max_mdp.vcf.gz --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/6_maxmissing0.8.vcf.gz
```

Need to bgzip to index the file
```
bcftools view ./qf4/2_quality_filtering/6_maxmissing0.8.vcf.gz -Oz -o ./qf4/2_quality_filtering/7_bgzip.vcf.gz
```

Create a bcftools index file for input into GATK
```
bcftools index ./qf4/2_quality_filtering/7_bgzip.vcf.gz
```

Create a GATK index file
```
gatk IndexFeatureFile -I ./qf4/2_quality_filtering/7_bgzip.vcf.gz --output ./qf4/2_quality_filtering/7_bgzip.vcf.gz.tbi
```
GATK standard quality filtering that only changes the filter column e.g. away from PASS 
```
gatk VariantFiltration \
-V ./qf4/2_quality_filtering/7_bgzip.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--verbosity ERROR \
-O ./qf4/2_quality_filtering/7b_gatk.vcf.gz
```

Check vcf file has chromosome names that are the same as the reference
```
zmore 5b_max_mdp.vcf.gz
```

Rename the chromosomes, where rename_chr.txt has col1: current name, col2:new name
```
bcftools annotate ./qf4/2_quality_filtering/7b_gatk.vcf.gz --rename-chrs ./qf4/2_quality_filtering/rename_chr.txt --output-type z --output ./qf4/2_quality_filtering/8_renamed.vcf.gz
```

Sort the file in chromosome and position order
```
bcftools sort ./qf4/2_quality_filtering/8_renamed.vcf.gz --output-type z --output ./qf4/2_quality_filtering/9_sorted.vcf.gz --temp-dir /g/data/ht96/Mel_UQ/qf4/2_quality_filtering/
```
Remove indels
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/9_sorted.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/10_SNPs_only.vcf.gz
```
		
Filtering for variants passing the filter
```
bcftools view ./qf4/2_quality_filtering/10_SNPs_only.vcf.gz -f.,PASS --output-type v --output ./qf4/2_quality_filtering/11_QF1.vcf
```

### QF2 (QF1 + MAF)
Minor allele frequency of 0.01
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/11_QF1.vcf.gz --maf 0.01 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/12_MAF.vcf.gz
```

### QF3 (QF1 + MAF + HWE)
Hardy Weinberg Equilibrium (HWE) 1e-6
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/12_QF2.vcf.gz --hwe 1e-6 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/13_QF3.vcf.gz
```

## 2. Population genetics analyses
### Thin and create PCA of QF2

Create list of thinned variants (thin.in)
```
./plink --vcf ./qf4/2_quality_filtering/12_QF2.vcf.gz \
--double-id --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex \
--indep-pairwise 50 10 0.1 \
--out ./qf4/3_popgen/4_thin_PCA/QF2_thin
```

Prune (using extract) and create pca
```
./plink --vcf ./qf4/2_quality_filtering/12_QF2.vcf.gz \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract ./qf4/3_popgen/4_thin_PCA/QF2_thin.prune.in --pca --out ./qf4/3_popgen/4_thin_PCA/QF2_PCA
```

### LD decay = poplddecay
Run PopLDdecay using QF2
```
./qf4/3_popgen/2_LDdecay/PopLDdecay/bin/PopLDdecay -InVCF ./qf4/2_quality_filtering/12_QF2.vcf -OutStat ./qf4/3_popgen/2_LDdecay/QF1_LDdecay 
```

### Nucleotide diversity (pi)
Output pi for 100kb windows on all sites
```
vcftools/bin/vcftools --vcf ./qf4/2_quality_filtering/11_QF1.vcf --window-pi 100000 --out ./qf4/3_popgen/5_AF_pi/100kb_pi_QF1
```
