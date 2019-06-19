##########Pré e pós imputação (Editado pela última vez em 04_06_19)
#PLINK v1.90b6.9 64-bit (4 Mar 2019)  
# por Jessica Mauer

############PARA SANGER#############
plink --bfile SCZ_QC_ibd --recode-vcf --keep-allele-order --out SCZ_paraimputar.vcf
bcftools +fixref SCZ_paraimputar.vcf -- -f human_g1k_v37.fasta
bcftools +fixref SCZ_paraimputar.vcf -Ov -o SCZ_paraimputar_strandflip.vcf -- -d -f human_g1k_v37.fasta -m flip
vcf-sort SCZ_paraimputar_strandflip.vcf | bgzip -c > vcfSCZ_fimSanger.vcf.gz

#########################
#### POST IMPUTATION ####
#########################

#######PARA SANGER#######
vi filterCHR1-8
i
bcftools filter 1.vcf.gz -e'MAF<0.01' -Oz -o 1filter.vcf.gz
bcftools filter 2.vcf.gz -e'MAF<0.01' -Oz -o 2filter.vcf.gz
bcftools filter 3.vcf.gz -e'MAF<0.01' -Oz -o 3filter.vcf.gz
bcftools filter 4.vcf.gz -e'MAF<0.01' -Oz -o 4filter.vcf.gz
bcftools filter 5.vcf.gz -e'MAF<0.01' -Oz -o 5filter.vcf.gz
bcftools filter 6.vcf.gz -e'MAF<0.01' -Oz -o 6filter.vcf.gz
bcftools filter 7.vcf.gz -e'MAF<0.01' -Oz -o 7filter.vcf.gz
bcftools filter 8.vcf.gz -e'MAF<0.01' -Oz -o 8filter.vcf.gz
ESC
:x

chmod 755 ./filterCHR1-8
nohup ./filterCHR1-8 &

vi filterCHR9-22
i
bcftools filter 9.vcf.gz -e'MAF<0.01' -Oz -o 9filter.vcf.gz
bcftools filter 10.vcf.gz -e'MAF<0.01' -Oz -o 10filter.vcf.gz
bcftools filter 11.vcf.gz -e'MAF<0.01' -Oz -o 11filter.vcf.gz
bcftools filter 12.vcf.gz -e'MAF<0.01' -Oz -o 12filter.vcf.gz
bcftools filter 13.vcf.gz -e'MAF<0.01' -Oz -o 13filter.vcf.gz
bcftools filter 14.vcf.gz -e'MAF<0.01' -Oz -o 14filter.vcf.gz
bcftools filter 15.vcf.gz -e'MAF<0.01' -Oz -o 15filter.vcf.gz
bcftools filter 16.vcf.gz -e'MAF<0.01' -Oz -o 16filter.vcf.gz
bcftools filter 17.vcf.gz -e'MAF<0.01' -Oz -o 17filter.vcf.gz
bcftools filter 18.vcf.gz -e'MAF<0.01' -Oz -o 18filter.vcf.gz
bcftools filter 19.vcf.gz -e'MAF<0.01' -Oz -o 19filter.vcf.gz
bcftools filter 20.vcf.gz -e'MAF<0.01' -Oz -o 20filter.vcf.gz
bcftools filter 21.vcf.gz -e'MAF<0.01' -Oz -o 21filter.vcf.gz
bcftools filter 22.vcf.gz -e'MAF<0.01' -Oz -o 22filter.vcf.gz
ESC
:x

chmod 755 ./filterCHR9-22
nohup ./filterCHR9-22 &

#juntando arquivos
bcftools concat -Oz 1filter.vcf.gz 2filter.vcf.gz 3filter.vcf.gz 4filter.vcf.gz 5filter.vcf.gz 6filter.vcf.gz 7filter.vcf.gz 8filter.vcf.gz 9filter.vcf.gz 10filter.vcf.gz 11filter.vcf.gz 12filter.vcf.gz 13filter.vcf.gz 14filter.vcf.gz 15filter.vcf.gz 16filter.vcf.gz 17filter.vcf.gz 18filter.vcf.gz 19filter.vcf.gz 20filter.vcf.gz 21filter.vcf.gz 22filter.vcf.gz -o scz_imp_bruta.vcf.gz

#Converter para .bed .bim .fam, excluindo variantes multialélicas
plink --vcf scz_imp_bruta.vcf.gz --biallelic-only strict --keep-allele-order --make-bed --out scz_imp

#Remover indels 
plink --bfile scz_imp --snps-only just-acgt --make-bed --out scz_imp_semindels

#Remover SNPs duplicadas
plink2 --bfile scz_imp_semindels --rm-dups exclude-all list --make-bed --out scz_imp_semduplicatas

#Fazendo arquivo updatesex:
awk '{print $1, $2, $5}' SCZ_QC_ibd.fam > upd_sex_scz.txt

#Fazendo arquivo fenotipo
awk '{print $1, $2, $6}' SCZ_QC_ibd.fam > upd_pheno_scz.txt

#Atualizando fenotipo e sexo
plink --bfile scz_imp_semduplicatas --update-sex upd_sex_scz.txt --pheno upd_pheno_scz.txt --make-bed --out scz_imputado_sex_pheno

#QC
plink --bfile scz_imputado_sex_pheno --maf 0.01 --geno 0.05 --hwe 0.001 --make-bed --out scz_imputado_QC




