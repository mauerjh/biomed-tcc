plink2 --bfile INPD_todos_qc --maf 0.01 --geno 0.1 --mind 0.1 --hwe 0.000001 --recode-vcf --out INPD_todos
#Separar os arquivos por cromossomos (chr1-chr22)
for i in `seq 1 22`; do plink2 --vcf INPD_todos.vcf --chr ${i} --recode-vcf --out chr${i}_INPD; done

for i in `seq 1 22`; do vcf-sort chr${i}_INPD.vcf | bgzip -c > chr${i}_INPD_sort1.vcf.gz; done

for i in '1 22'; do bcftools index chr${i}_INPD_sort1.vcf.gz ; done

for i in 'seq 1 22'; do bcftools norm -cx chr${i}_INPD_sort1.vcf.gz -f human_g1k_v37.fasta -Oz chr${i}_INPD_TOPMED_ref.vcf.gz; done

#Create a sorted *.vcf.gz file using VCFtools and tabix (including bgzip):
for i in `seq 1 22`; do vcf-sort chr${i}_INPD_TOPMED_ref.vcf| bgzip -c > chr${i}_INPD_fimTOPMED.vcf.gz; done

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
bcftools concat -Oz 1filterinfo.vcf.gz 2filterinfo.vcf.gz 3filterinfo.vcf.gz 4filterinfo.vcf.gz 5filterinfo.vcf.gz 6filterinfo.vcf.gz 7filterinfo.vcf.gz 8filterinfo.vcf.gz 9filterinfo.vcf.gz 10filterinfo.vcf.gz 11filterinfo.vcf.gz 12filterinfo.vcf.gz 13filterinfo.vcf.gz 14filterinfo.vcf.gz 15filterinfo.vcf.gz 16filterinfo.vcf.gz 17filterinfo.vcf.gz 18filterinfo.vcf.gz 19filterinfo.vcf.gz 20filterinfo.vcf.gz 21filterinfo.vcf.gz 22filterinfo.vcf.gz -o imp_bruta.vcf.gz

#### Merge results GWAS with imputed information file and filter by R2 >= 6 (the results have been filtered by MAF before)
# Extract columns of interest in GWAS results: CHR POS SNPID Allele1 Allele2 AF_Allele2 N BETA SE Tstatp.value
zcat AD.chr19.SAIGE.txt.gz | awk '{print
$1,$2,$3,$4,$5,$7,$8,$9,$10,$11,$12}' > temp_AD.chr19.results.txt
# Extract columns of interest in imputation information: SNP Rsq
zcat ../saige/chr19.info.gz | awk '{print $1,$7}' > temp_chr19.info.txt
# Merge information from the GWAS and R2
awk 'NR==FNR{a[$1]=$2; next} $3 in a{print $0,a[$3]}' temp_chr19.info.txt temp_AD.chr19.results.txt > temp
# Filer by R2 >= 0.6 and put the header
echo 'CHR POS SNPID Allele1 Allele2 AF_Allele2 N BETA SE Tstat p.value
Rsq' > AD.chr19.results.QC.txt
awk '{if ($12>=0.6) print $0}' temp >> AD.chr19.results.QC.txt
#### Plot the results: Manhattan plot and regional plot
awk '{if (NR>1) print $1,$2,$11}' temp_AD.chr19.results.QC.txt >
temp_plot.AD.chr19.txt
sort -k11 -g temp_AD.chr19.results.QC.txt | head # lowest p-value at
19:45422946 , p = 9.3515193823689e-16
awk '{if (NR==1) print $1,$2,$3,$11; else if ($2>= 44922946 && $2<=
45922946) print $1,$2,"chr"$3,$11}' AD.chr19.results.QC.txt >
temp_ld.AD.chr19.txt


#Converter para plink2 format (pgen, pvar, psam)
plink2 --vcf imp_bruta.vcf.gz --make-pgen --out imputacao_bruta

#Remover SNPs duplicadas
plink2 --pfile imputacao_bruta --rm-dup exclude-all list --make-pgen --out imp_semduplicatas

#Remover indels 
plink2 --pfile imp_semduplicatas --snps-only just-acgt --make-pgen --out imp_semindels

#Fazendo arquivo updatesex:
awk '{print $2, $5}' arquivo_antes_de_imputar.fam > upd_sex.txt
gedit upd_sex.txt
#colocar cabeçalho: '#IID SEX', salvar
#Fazendo arquivo fenotipo
awk '{print $2, $6}' arquivo_antes_de_imputar.fam > upd_pheno.txt

#Atualizando fenotipo e sexo
plink2 --pfile imp_semindels --update-sex upd_sex.txt --pheno upd_pheno.txt --make-pgen --out imputado_sex_pheno

#QC
plink2 --pfile imputado_sex_pheno --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000001 --make-pgen --out imputado_QC
