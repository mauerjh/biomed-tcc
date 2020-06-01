############PARA SANGER#############

../../plink --bfile pep_qc_final --maf 0.01 --geno 0.1 --mind 0.1 --hwe 0.000001 --recode-vcf --out PEP
vcf-sort PEP.vcf | bgzip -c > vcfPEP_sort.vcf.gz
bcftools norm -cx vcfPEP_sort.vcf.gz -f human_g1k_v37.fasta -o vcfPEP_ref.vcf

for i in `seq 1 22`; do ../../plink --vcf vcfPEP_ref.vcf --chr ${i} --recode-vcf --out chr${i}_PEP; done
for i in `seq 1 22`; do vcf-sort hgdp_chr${i}.vcf | bgzip -c > hgdp_chr${i}.vcf.gz; done

############PARA MICHIGAN#############

plink --bfile INPD_todos_qc --recode-vcf --out INPD_todos
#Separar os arquivos por cromossomos (chr1-chr22)
for i in `seq 1 22`; do plink2 --vcf INPD_todos.vcf --chr ${i} --recode-vcf --out chr${i}_INPD; done

for i in `seq 1 22`; do vcf-sort chr${i}_INPD.vcf | bgzip -c > chr${i}_INPD_sort1.vcf.gz; done

for i in '1 22'; do bcftools index chr${i}_INPD_sort1.vcf.gz ; done

for i in 'seq 1 22'; do bcftools norm -cx chr${i}_INPD_index.vcf.gz -f human_g1k_v37.fasta -Oz chr${i}_INPD_TOPMED_ref.vcf.gz; done

#Create a sorted *.vcf.gz file using VCFtools and tabix (including bgzip):
for i in `seq 1 22`; do vcf-sort chr${i}_INPD_TOPMED_ref.vcf| bgzip -c > chr${i}_INPD_fimTOPMED.vcf.gz; done

ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/human_g1k_v37.fasta.fai
