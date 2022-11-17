# MangoWGS
225 Mangifera indica, raw vcf file with non-variants

1. SNP filtering from raw vcf
2. Population genetics analyses
3. Identifying deleterious mutations using SIFT

## SNP filtering

### minGQ, Minimum genotype quality of 20 

```
vcftools --gzvcf ./qf4/1_rawVCF/4f90968ba585820a0142bd5dc0d5562d_raw_variants_1.vcf.gz --minGQ 20 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/1_GQ.recode.vcf.gz
```

### minDP, Minimum depth per sample of 5
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/1_GQ.recode.vcf.gz --minDP 5 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/2_mindp.vcf.gz
```

### Missing data per site: initially a relaxed threshold (i.e. 50%)
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/2_mindp.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/3_max-missing.vcf.gz
```

### Create a list of missing data per individual:
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/3_max-missing.vcf.gz --missing-indv --out ./qf4/2_quality_filtering/4_
```
No individuals were removed as all have more than 85% of data
  
## Maximum mean depth for removal of paralogues
### Create a list of mean depth per site:
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/3_max-missing.vcf.gz --site-mean-depth --out ./qf4/2_quality_filtering/5_mean_depth
```
### Filter for maximum mean depth of 50 for all sites
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/3_max-missing.vcf.gz --max-meanDP 50 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/5b_max_mdp.vcf.gz
```	


### More stringent missing data per locus, filter for overall 20% missing data:
```
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/5b_max_mdp.vcf.gz --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/6_maxmissing0.8.vcf.gz
```

