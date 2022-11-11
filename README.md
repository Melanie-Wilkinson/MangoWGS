# MangoWGS
## 225 Mangifera indica

### Raw vcf file with non-variants:
./Q4815/AS17000_joint_raw_indica_20221027/4f90968ba585820a0142bd5dc0d5562d_raw_variants_1.vcf.gz

### minGQ, Minimum genotype quality of 20
vcftools --gzvcf ./qf4/1_rawVCF/4f90968ba585820a0142bd5dc0d5562d_raw_variants_1.vcf.gz --minGQ 20 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/1_GQ.recode.vcf.gz
    ### Memory used 94.69GB, wall time used: 25:54:40hrs

### minDP, Minimum depth per sample of 5
vcftools/bin/vcftools --gzvcf ./qf4/2_quality_filtering/1_GQ.recode.vcf.gz --minDP 5 --recode --recode-INFO-all --stdout | gzip -c > ./qf4/2_quality_filtering/2_mindp.vcf.gz
	    # Memory used <95GB, wall time used: 31:18:25