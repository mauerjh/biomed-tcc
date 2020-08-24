#################
# TOPMED R2 IMPUTATION (MINIMAC4)
# INPUT arquivos pós qc 
#################

../plink --bfile INPD2020_sempaisfam_semXY --maf 0.01 --geno 0.1 --mind 0.1 --hwe 0.000001 --recode vcf --keep-allele-order --out INPD2020_qc
for i in `seq 1 22`; do ../plink --vcf INPD2020_qc.vcf --chr ${i} --recode vcf --keep-allele-order --out chr${i}_INPD2020; done
for i in `seq 1 22`; do vcf-sort chr${i}_INPD2020.vcf | bgzip -c > chr${i}_INPD2020_sort1.vcf.gz; done
for i in {1..22}; do echo bcftools index chr${i}_INPD2020_sort1.vcf.gz ; done >> indexar
chmod 755 indexar
./indexar
for i in {1..22}; do echo bcftools norm -cx chr${i}_INPD2020_sort1.vcf.gz -f human_g1k_v37.fasta -o chr${i}_INPD2020_TOPMED_ref.vcf; done >> normalizar
chmod 755 normalizar
nohup ./normalizar &
for i in `seq 1 22`; do vcf-sort chr${i}_INPD2020_TOPMED_ref.vcf| bgzip -c > chr${i}_INPD2020_fimTOPMED.vcf.gz; done
 
ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai

#########################
#### POST IMPUTATION ####
#########################
#Filtrar variantes por MAF

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
bcftools concat -Oz 1filter.vcf.gz 2filter.vcf.gz 3filter.vcf.gz 4filter.vcf.gz 5filter.vcf.gz 6filter.vcf.gz 7filter.vcf.gz 8filter.vcf.gz 9filter.vcf.gz 10filter.vcf.gz 11filter.vcf.gz 12filter.vcf.gz 13filter.vcf.gz 14filter.vcf.gz 15filter.vcf.gz 16filter.vcf.gz 17filter.vcf.gz 18filter.vcf.gz 19filter.vcf.gz 20filter.vcf.gz 21filter.vcf.gz 22filter.vcf.gz -o imp_bruta.vcf.gz

#Converter para plink2 format (pgen, pvar, psam)
plink2 --vcf imp_bruta.vcf.gz --make-pgen --out imputacao_bruta

#Remover SNPs duplicadas
plink2 --pfile imputacao_bruta --rm-dup exclude-all list --make-pgen --out imp_semduplicatas

#Remover indels 
plink2 --pfile imp_semduplicatas --snps-only just-acgt --make-bed --out imp_semindels

#Fazendo arquivo updatesex:
awk '{print 0, $2, $5}' arquivo_antes_de_imputar.fam > upd_sex.txt
paste upd_sex.txt imp_semindels.fam | awk '{print $1, $2, $6,$7.$3,$9}' > imputado_sex.fam
cp imp_semindels.bed imputado_sex.bed
cp imp_semindels.bim imputado_sex.bim
r
#Fazendo arquivo fenotipo
awk '{print $2, $6}' arquivo_antes_de_imputar.fam > upd_pheno.txt

#QC
plink2 --bfile imputado_sex --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000001 --make-bed --out imputado_QC
