##########################
# Script para plotar pcs em cima do 1000genomes, colorindo por coorte
# adaptado de
# Editado pela ultima vez em 05/09/2019
# Jessica Mauer
##########################

# Baixar dados 1K genomes

wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.bed.gz"
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.bim.gz"
wget "http://genepi.qimr.edu.au/staff/sarahMe/enigma/MDS/HM3_b37.fam.gz"  

# QC e merge do dataset de interesse com 1k genomes

export datafileraw=merge_semasiaticos_final

./plink --bfile $datafileraw --hwe 1e-6 --geno 0.03 --maf 0.01 --make-bed --out ${datafileraw}_filtered
 
gunzip HM3_b37*.gz

export datafile=${datafileraw}_filtered

awk '{print $2}' HM3_b37.bim > HM3_b37.snplist.txt

./plink --bfile $datafile --extract HM3_b37.snplist.txt --make-bed --out local

awk '{if (($5=="T" && $6=="A")||($5=="A" && $6=="T")||($5=="C" && $6=="G")||($5=="G" && $6=="C")) print $2, "ambig"; else print $2;}' $datafile.bim | grep -v ambig > local.snplist.txt

../plink --bfile HM3_b37 --extract local.snplist.txt --make-bed --out external

../plink --bfile local --bmerge external.bed external.bim external.fam --make-bed --out HM3_b37merge

../plink --bfile local --flip HM3_b37merge-merge.missnp --make-bed --out flipped

../plink --bfile flipped --bmerge external.bed external.bim external.fam --make-bed --out HM3_b37merge

../plink --bfile flipped --exclude HM3_b37merge-merge.missnp --make-bed --out flip_excluded

../plink --bfile flip_excluded --bmerge external.bed external.bim external.fam --make-bed --out HM3_b37merge

# Cálculo PCs

../plink --bfile HM3_b37merge_semGIH --cluster --mind .05 --mds-plot 4 --extract local.snplist.txt --out HM3_b37mds

awk 'BEGIN{OFS=","};{print $1,$2,$3,$4,$5,$6,$7}' >> HM3_b37mds2R.mds.csv HM3_b37mds.mds 

###### Plot #########

# carregar e arrumar dados

R

library(calibrate)

mds.cluster = read.csv("HM3_b37mds2R.mds.csv", header = T)

colors=rep(vector(),length(mds.cluster$C1))

colors[which(mds.cluster$FID == "CEU")] <- "lightblue"

colors[which(mds.cluster$FID == "CHB")] <- "brown"

colors[which(mds.cluster$FID == "YRI")] <- "yellow"

colors[which(mds.cluster$FID == "TSI")] <- "green4"

colors[which(mds.cluster$FID == "JPT")] <- "purple"

colors[which(mds.cluster$FID == "CHD")] <- "orange"

colors[which(mds.cluster$FID == "MEX")] <- "grey50"

colors[which(mds.cluster$FID == "ASW")] <- "darkolivegreen"

colors[which(mds.cluster$FID == "LWK")] <- "magenta"

colors[which(mds.cluster$FID == "MKK")] <- "darkblue"

# Renomear FIDs de acordo com a coorte, mudar numeros das linhas de acordo com o seu dataset


mds.cluster[2753:2926, "FID"] <- "SCZ" #T
mds.cluster[2420:2183, "FID"] <- "SCZ" #C
mds.cluster[2694:2752, "FID"] <- "PEP" #P e PW
mds.cluster[2574:2643, "FID"] <- "PEP" #CP
mds.cluster[2421:2533, "FID"] <- "PEP" #C
mds.cluster[2534:2573, "FID"] <- "TEPT" #CNP
mds.cluster[2644:2693, "FID"] <- "TEPT" #NP

# Obs, eu salvei um arquivo com os FIDs correspondentes as coortes
# write.csv(mds.cluster, file = "HM3_b37mds2R_coortes.mds.csv", quote = F, row.names = F, col.names = T)
# Colorir coortes - o inpd vai sobrar em vermelho

colors[which(mds.cluster$FID == "PEP")] <- "green"

colors[which(mds.cluster$FID == 0)] <- "green"

colors[which(mds.cluster$FID == "SCZ")] <- "blue"

colors[which(mds.cluster$FID == "TEPT")] <- "cyan"

# Plot de fato e salvar em pdf

pdf(file="mdsplottestec1xc2.pdf", width=7, height=7)
plot(rev(mds.cluster$C2), rev(mds.cluster$C1), col=rev(colors), ylab="Dimension 1", xlab="Dimension 2", pch=20, cex = 0.5)
legend("topright", c("INPD", "PEP", "SCZ", "TEPT", "CEU","CHB","YRI","TSI","JPT","CHD","MEX","GIH","ASW","LWK","MKK"), fill=c("red", "green", "blue", "cyan", "lightblue", "brown", "yellow", "green4","purple","orange","grey50","black","darkolivegreen","magenta","darkblue"))
dev.off()

pdf(file="mdsplotc1xc3.pdf", width=7, height=7)
plot(rev(mds.cluster$C3), rev(mds.cluster$C1), col=rev(colors), ylab="Dimension 1", xlab="Dimension 3", pch=20, cex = 0.5)
legend("topright", c("INPD", "PEP", "SCZ", "TEPT", "CEU","CHB","YRI","TSI","JPT","CHD","MEX","GIH","ASW","LWK","MKK"), fill=c("red", "springgreen4", "blue", "cyan", "lightblue", "brown", "yellow", "green","purple","orange","grey50","black","darkolivegreen","magenta","darkblue"))
dev.off()

pdf(file="mdsplotc2xc3.pdf", width=7, height=7)
plot(rev(mds.cluster$C3), rev(mds.cluster$C2), col=rev(colors), ylab="Dimension 2", xlab="Dimension 3", pch=20, cex = 0.5)
legend("topright", c("INPD", "PEP", "SCZ", "TEPT", "CEU","CHB","YRI","TSI","JPT","CHD","MEX","GIH","ASW","LWK","MKK"), fill=c("red", "springgreen4", "blue", "cyan", "lightblue", "brown", "yellow", "green","purple","orange","grey50","black","darkolivegreen","magenta","darkblue"))
dev.off()

# plot 3d

library(calibrate)
library(rgl)
plot3d(rev(mds.cluster$C1), rev(mds.cluster$C2), rev(mds.cluster$C3),col=rev(colors))


