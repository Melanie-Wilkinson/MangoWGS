# MangoWGS
225 Mangifera indica

1. Alignments and variant calling

2. Quality filtering

3. Centromere estimation

4. Local Principal Component Analysis
   
5. Population genetics analyses	 <br>
	5a. Homozygosity <br>
	5b. Genetic differentiation (FST) <br>
	5c. Nucleotide diversity (pi) <br>
	5d. Linkage disequilibrium (LD) heatmap <br> 

6. Genome-wide association study (GWAS)

8. Identifying deleterious mutations using SIFT


## 1. Alignments and variant calling
   


## 2. SNP filtering
### QF1
minGQ, Minimum genotype quality of 20 using ```VCFtools v0.1.17``` [(Danecek et al. 2011)](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).

```
vcftools --gzvcf 4f90968ba585820a0142bd5dc0d5562d_raw_variants_1.vcf.gz --minGQ 20 --recode --recode-INFO-all --stdout | gzip -c > 1_GQ.recode.vcf.gz
```

minDP, Minimum depth per sample of 5 using ```VCFtools v0.1.17```.
```
vcftools --gzvcf 1_GQ.recode.vcf.gz --minDP 5 --recode --recode-INFO-all --stdout | gzip -c > 2_mindp.vcf.gz
```

Missing data per site: initially a relaxed threshold (i.e. 50%) using ```VCFtools v0.1.17```.
```
vcftools --gzvcf 2_mindp.vcf.gz --max-missing 0.5 --recode --recode-INFO-all --stdout | gzip -c > 3_max-missing.vcf.gz
```

Create a list of missing data per individual using ```VCFtools v0.1.17```:
```
vcftools --gzvcf 3_max-missing.vcf.gz --missing-indv --out 4_
```
No individuals were removed as all have more than 85% of data
  
Maximum mean depth for removal of paralogues
1. Create a list of mean depth per site using ```VCFtools v0.1.17```:
```
vcftools --gzvcf 3_max-missing.vcf.gz --site-mean-depth --out 5_mean_depth
```
2. Filter for maximum mean depth of 50 for all sites using ```VCFtools v0.1.17```.
```
vcftools --gzvcf 3_max-missing.vcf.gz --max-meanDP 50 --recode --recode-INFO-all --stdout | gzip -c > 5b_max_mdp.vcf.gz
```	


More stringent missing data per locus, filter for overall 20% missing data using ```VCFtools v0.1.17```:
```
vcftools --gzvcf 5b_max_mdp.vcf.gz --max-missing 0.8 --recode --recode-INFO-all --stdout | gzip -c > 6_maxmissing0.8.vcf.gz
```

Bgzip the index file using ```bcftools v1.12```  [(Danecek et al. 2021)](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).
```
bcftools view 6_maxmissing0.8.vcf.gz -Oz -o 7_bgzip.vcf.gz 
```

Create a bcftools index file for input into GATK using ```bcftools v1.12```.
```
bcftools index 7_bgzip.vcf.gz
```

Create a GATK index file using ```GATK v4.2.5``` [(Van der Auwera & O'Connor 2020)](https://www.oreilly.com/library/view/genomics-in-the/978149197518).
```
gatk IndexFeatureFile -I 7_bgzip.vcf.gz --output 7_bgzip.vcf.gz.tbi 
```
```GATK v4.2.5``` standard quality filtering that only changes the filter column e.g. away from PASS 
```
gatk VariantFiltration \
-V 7_bgzip.vcf.gz \
-filter "QD < 2.0" --filter-name "QD2" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--verbosity ERROR \
-O 7b_gatk.vcf.gz
```

Check vcf file has chromosome names that are the same as the reference
```
zmore 5b_max_mdp.vcf.gz
```

Rename the chromosomes, where rename_chr.txt has col1: current name, col2:new name using ```bcftools v1.12```.
```
bcftools annotate 7b_gatk.vcf.gz --rename-chrs rename_chr.txt --output-type z --output 8_renamed.vcf.gz
```

Sort the file in chromosome and position order using ```bcftools v1.12```.
```
bcftools sort 8_renamed.vcf.gz --output-type z --output 9_sorted.vcf.gz
```
Remove indels using ```VCFtools v0.1.17```.
```
vcftools --gzvcf 9_sorted.vcf.gz --remove-indels --recode --recode-INFO-all --stdout | gzip -c > 10_SNPs_only.vcf.gz
```
		
Filtering for variants passing the filter using ```bcftools v1.12```.
```
bcftools view 10_SNPs_only.vcf.gz -f.,PASS --output-type v --output 11_QF1.vcf
```

### QF2 (QF1 + MAF)
Remove sites with a minor allele frequency of < 0.01 using ```VCFtools v0.1.17```.
```
vcftools --gzvcf 11_QF1.vcf.gz --maf 0.01 --recode --recode-INFO-all --stdout | gzip -c > 12_MAF.vcf.gz
```

## 3. Centromere estimation 
Estimate centromere location in _Mangifera indica_ cv. 'Alphonso' (version: CATAS_Mindica_2.1) using  ```RepeatOBserverV1``` [(Elphinstone et al. 2023)](https://www.biorxiv.org/content/10.1101/2023.12.30.573697v1) <br>
Alphonso genome info: <br>
 https://www.ncbi.nlm.nih.gov/datasets/taxonomy/29780/ <br>
 https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011075055.1/

Download the Setup_Run_Repeats.sh from github
```
wget https://raw.githubusercontent.com/celphin/RepeatOBserverV1/main/Setup_Run_Repeats.sh
chmod +x Setup_Run_Repeats.sh
dos2unix Setup_Run_Repeats.sh
```

To install the R package "RepeatOBserverV1", you will first need to install the package devtools in your version of R.
```
install.packages("devtools")
library(devtools)
install_github("celphin/RepeatOBserverV1")
library(RepeatOBserverV1)
```

------------------------------
go to FTP directory for GenBank assembly (version: CATAS_Mindica_2.1)
```
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/011/075/055/GCF_011075055.1_CATAS_Mindica_2.1/GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna.gz

gunzip GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna.gz

mv GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna Mango.fasta
```
---------------------------------
Estimate centromeres
```
module load StdEnv/2020
module load seqkit/2.3.1
module load emboss/6.6.0
module load r/4.1.2

srun Setup_Run_Repeats.sh -i Mango -f Mango.fasta -h H0 -c 20 -m 191000M -g FALSE
```

## 4. Local Principal Component Analysis 

Remove invariant sites using ```GATK v4.2.5``` SelectVariants. 
```
gatk SelectVariants -R GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna -V 11_QF1.vcf.gz --exclude-non-variants --select-type-to-include SNP -O QF1_variant_only.vcf 
```
Remove contigs that are not aligned to a chromosome using ```bcftools v1.12```. 
```
bcftools view ${DIR}/11_QF1_variant_only.vcf.gz -R chrom_list.tsv -o QF1_variant_only.vcf 
```
Remove sites with greater than 95% missing data using ```VCFtools v0.1.17```.
```
vcftools --gzvcf 11_QF1_variant_only.vcf.gz --max-missing 0.95 --recode --recode-INFO-all --out QF1_variant_only_95miss
```
Convert SNP dataset to BCF format using ```bcftools v1.12```.
```
bcftools convert -O b QF1_variant_only_95miss.vcf > QF1_variant_only_95miss.bcf
bcftools index QF1_variant_only_95miss.bcf
```
Perform local PCA with the R package ```lostruct v0.0.0.9000``` [(Li & Ralph 2019)](https://academic.oup.com/genetics/article/211/1/289/5931130?login=true) on windows of 1000SNPs.
Perform PCA of outlier regions identified through local PCA analysis. 
These analyses was performed using code from and as detailed in [(Huang et al. 2020)](https://onlinelibrary-wiley-com.ezproxy.library.uq.edu.au/doi/10.1111/mec.15428). 

## 5. Inversions - population genetics analyses
### 5a. Homozygosity/heterozygosity 
Output heterozygosity in 10kb windows in miinv3.0 for every individual
```
vcftools --gzvcf 12_QF2.vcf.gz --het --chr NC_058139.1 --out 1_chr3 --from-bp 14606714 --to-bp 17874417 
```

### 5b. Genetic differentiation (Fst)
Output Fst in 10kb windows between cluster 0 and 2 from the local PCA for miinv3.0
```
vcftools --gzvcf 12_QF2.vcf.gz --weir-fst-pop 1_chr3_0.txt --weir-fst-pop 1_chr3_2.txt --fst-window-size 10000 --chr NC_058139.1 --out 1_chr3_FST
```

### 5c. Nucleotide diversity (pi)
Output pi in 10kb windows for each cluster from the local PCA for miinv3.0 using pixy
```
pixy --stats pi --vcf 11_QF1.vcf.gz --window_size 10000 --populations 1_chr3_pop.txt --chromosomes 'NC_058139.1' --output_prefix '1_chr3'
```

### 5d. Linkage disequilibrium (LD) heatmap
We used the script emerald2windowldcounts.pl from https://github.com/owensgl/reformat <br>
LD for every chromosome seperately using a loop
```
for chr in "${chr[@]}"; do
bcftools view -r ${chr} 12_QF2bgzip.vcf.gz | vcftools --vcf - --maf 0.05 --thin 1000 -c --geno-r2 --max-missing-count 0 | perl emerald2windowldcounts.pl > ${chr}.ld.txt

done
```

## 6. Genome-wide association study (GWAS) on continuous phenotype

Create list of thinned variants (thin.in) for chr1-20, where '--indep-pairwise' is consider a window of 50kb, calculate LD between each pair of SNPs in the window, remove one of a pair of SNPs if the LD is greater than 0.1, shift the window 10 SNPs forward and repeat the procedure using ```plink v1.9``` [(Purcell et al. 2009)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1950838/) 

```
plink --vcf 12_QF2.vcf.gz \
--double-id --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex \
--indep-pairwise 50 10 0.1 \
--chr NC_058137.1, NC_058138.1, NC_058139.1, NC_058140.1, NC_058141.1, NC_058142.1, NC_058143.1, NC_058144.1, NC_058145.1, NC_058146.1, NC_058147.1, NC_058148.1, NC_058149.1, NC_058150.1, NC_058151.1, NC_058152.1, NC_058153.1, NC_058154.1, NC_058155.1, NC_058156.1 \
--out QF2_thin
```

Prune (using extract) and create PCA (eigenvec and eigenval used as covariants in GWAS) using ```plink v1.9```.
```
plink --vcf 12_QF2.vcf.gz \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--extract QF2_thin.prune.in --pca --out QF2_pruned
```

Association test where phenotype file is ID1 ID2 phenotype_value and PC1 - PC6 are covariants using ```plink v1.9```.
```
plink --vcf 12_QF2.vcf.gz --pheno TC_Age12_WRS_GWAS.txt \
--double-id --allow-extra-chr --set-missing-var-ids @:# --allow-no-sex \
--adjust qq-plot log10 --ci 0.95 --covar QF2_pruned.eigenvec --covar-number 1-6 --linear \
--extract QF2_thin.prune.in --out QF2_pruned
```

## 6. Identifying deleterious mutations
Build a genome database for mango using [Sorting Intolerant From Tolerant 4G create genomic database](https://github.com/pauline-ng/SIFT4G_Create_Genomic_DB)

Download Alphonso genome file: GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna.gz (last modified 2021-10-19 22:28)
```
cd chr-src
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Mangifera_indica/latest_assembly_versions/GCF_011075055.1_CATAS_Mindica_2.1/GCF_011075055.1_CATAS_Mindica_2.1_genomic.fna.gz'
```
Download compressed annotation file: GCF_011075055.1_CATAS_Mindica_2.1_genomic.gtf.gz 
```
cd gene-annotation-src
wget 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/plant/Mangifera_indica/latest_assembly_versions/GCF_011075055.1_CATAS_Mindica_2.1/GCF_011075055.1_CATAS_Mindica_2.1_genomic.gtf.gz'
```

Create a config file  
Use test_files/homo_sapiens-test.txt as a template and modify as below <br>

```
		GENETIC_CODE_TABLE=1
		GENETIC_CODE_TABLENAME=Standard

		PARENT_DIR={DIR}/mango_genome		
		ORG=mangifera_indica
		ORG_VERSION=assembly_CATAS_Mindica_2.1

		#Running SIFT 4G
		SIFT4G_PATH=/sift4g/bin/sift4g
		PROTEIN_DB=/SIFT_database/uniref90.fasta

		**Removed** - Mitochondria or cholorplast (plastid) annotation files and DBSNP_VCF_FILE 
```

Create database
```
perl make-SIFT-db-all.pl -config mango_config.txt
```	


Find deleterious scores for all sites using SIFT 4G [(Vaser et al. 2016)](https://www.nature.com/articles/nprot.2015.123) 
```
java -jar SIFT4G_Annotator.jar -c -i 11_QF1.vcf -d assembly_CATAS_Mindica_2.1 -r 1_225i_all
```

Output allelic frequency for specified sites using ```VCFtools v0.1.17```.
```
vcftools --vcf 11_QF1.vcf --freq --positions del_mutations.txt
```		

Output genotypes for specified sites using ```bcftools v1.12```.
```
bcftools query -H -f '%CHROM:%POS\t%REF\t%ALT[\t%GT]\n' 11_QF1.vcf.gz --regions NC_058141.1:7109828
```

