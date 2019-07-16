##########Pré e pós imputação Sanger (ultima atualizacao em 06_06_19)
#PLINK v2.00a2LM 64-bit Intel (27 May 2019) 
#bcftools Version: 1.4.1-1-gcdac3db (using htslib 1.4.1-1-g7771288)
#VCFtools (0.1.15)
# por Jessica Mauer

############PARA SANGER#############

# para arquivos saídos do qc feito em plink 1.9:
plink2 --bfile SCZ_QC_ibd --make-pgen --out SCZ_imput

## para arquivos saídos do qc feito em plink 2.0, começar aqui:
plink2 --pfile SCZ_imput --export vcf --out SCZ_imput_vcf

#stats:
bcftools +fixref SCZ_imput_vcf.vcf -- -f human_g1k_v37.fasta

#correçao alelo referencia com base no fasta b37/hg19
bcftools +fixref SCZ_imput_vcf.vcf -Ov -o SCZ_paraimputar_strandflip.vcf -- -d -f human_g1k_v37.fasta -m flip
vcf-sort SCZ_paraimputar_strandflip.vcf | bgzip -c > vcfSCZ_fimSanger.vcf.gz

#########################
#### POST IMPUTATION ####
#########################

#######PARA SANGER#######
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

#Filtrar variantes por qualidade da imputação (info score)

vi filterINFO1-8
i
bcftools filter 1filter.vcf.gz -e 'INFO<0.8' -Oz -o 1filterinfo.vcf.gz
bcftools filter 2filter.vcf.gz -e 'INFO<0.8' -Oz -o 2filterinfo.vcf.gz
bcftools filter 3filter.vcf.gz -e 'INFO<0.8' -Oz -o 3filterinfo.vcf.gz
bcftools filter 4filter.vcf.gz -e 'INFO<0.8' -Oz -o 4filterinfo.vcf.gz
bcftools filter 5filter.vcf.gz -e 'INFO<0.8' -Oz -o 5filterinfo.vcf.gz
bcftools filter 6filter.vcf.gz -e 'INFO<0.8' -Oz -o 6filterinfo.vcf.gz
bcftools filter 7filter.vcf.gz -e 'INFO<0.8' -Oz -o 7filterinfo.vcf.gz
bcftools filter 8filter.vcf.gz -e 'INFO<0.8' -Oz -o 8filterinfo.vcf.gz
ESC
:x

chmod 755 ./filterINFO1-8
nohup ./filterINFO1-8 &

vi filterINFO9-22
i
bcftools filter 9filter.vcf.gz -e 'INFO<0.8' -Oz -o 9filterinfo.vcf.gz
bcftools filter 10filter.vcf.gz -e 'INFO<0.8' -Oz -o 10filterinfo.vcf.gz
bcftools filter 11filter.vcf.gz -e 'INFO<0.8' -Oz -o 11filterinfo.vcf.gz
bcftools filter 12filter.vcf.gz -e 'INFO<0.8' -Oz -o 12filterinfo.vcf.gz
bcftools filter 13filter.vcf.gz -e 'INFO<0.8' -Oz -o 13filterinfo.vcf.gz
bcftools filter 14filter.vcf.gz -e 'INFO<0.8' -Oz -o 14filterinfo.vcf.gz
bcftools filter 15filter.vcf.gz -e 'INFO<0.8' -Oz -o 15filterinfo.vcf.gz
bcftools filter 16filter.vcf.gz -e 'INFO<0.8' -Oz -o 16filterinfo.vcf.gz
bcftools filter 17filter.vcf.gz -e 'INFO<0.8' -Oz -o 17filterinfo.vcf.gz
bcftools filter 18filter.vcf.gz -e 'INFO<0.8' -Oz -o 18filterinfo.vcf.gz
bcftools filter 19filter.vcf.gz -e 'INFO<0.8' -Oz -o 19filterinfo.vcf.gz
bcftools filter 20filter.vcf.gz -e 'INFO<0.8' -Oz -o 20filterinfo.vcf.gz
bcftools filter 21filter.vcf.gz -e 'INFO<0.8' -Oz -o 21filterinfo.vcf.gz
bcftools filter 22filter.vcf.gz -e 'INFO<0.8' -Oz -o 22filterinfo.vcf.gz
ESC
:x

chmod 755 ./filterINFO9-22
nohup ./filterINFO9-22 &

#juntando arquivos
bcftools concat -Oz 1filterinfo.vcf.gz 2filterinfo.vcf.gz 3filterinfo.vcf.gz 4filterinfo.vcf.gz 5filterinfo.vcf.gz 6filterinfo.vcf.gz 7filterinfo.vcf.gz 8filterinfo.vcf.gz 9filterinfo.vcf.gz 10filterinfo.vcf.gz 11filterinfo.vcf.gz 12filterinfo.vcf.gz 13filterinfo.vcf.gz 14filterinfo.vcf.gz 15filterinfo.vcf.gz 16filterinfo.vcf.gz 17filterinfo.vcf.gz 18filterinfo.vcf.gz 19filterinfo.vcf.gz 20filterinfo.vcf.gz 21filterinfo.vcf.gz 22filterinfo.vcf.gz -o imp_bruta.vcf.gz

#Converter para plink2 format (pgen, pvar, psam)
plink2 --vcf scz_imp_bruta.vcf.gz --make-pgen --out scz_imp

#Remover SNPs duplicadas
plink2 --bfile scz_imp_semindels --rm-dups exclude-all list --make-bed --out scz_imp_semduplicatas

#Remover indels 
plink2 --pfile scz_imp --snps-only just-acgt --make-pgen --out scz_imp_semindels

#Fazendo arquivo updatesex:
awk '{print $2, $5}' SCZ_imput.psam > upd_sex_scz.txt

#Fazendo arquivo fenotipo
awk '{print $2, $6}' SCZ_imput.psam > upd_pheno_scz.txt

#Atualizando fenotipo e sexo
plink2 --pfile scz_imp_semindels --update-sex upd_sex_scz.txt --pheno upd_pheno_scz.txt --make-pgen --out scz_imputado_sex_pheno

#QC
plink2 --pfile scz_imputado_sex_pheno --maf 0.01 --mind 0.1 --geno 0.1 --hwe 0.000001 --make-pgen --out scz_imputado_QC




