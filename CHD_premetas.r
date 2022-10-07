#Transcribed script from main CHD script for the generation of the premeta files going into the final analysis
#3/12/18: So far we have the following studies:
#UKBB - HAIL adjusted
#EPIC
#Cardiogram - 1KG
#Cardiogram - Exome
#Cardiogram - Metabo
#GerMIFs1-3,4-7
#deCODE
#Greek
#HUNT
#Partners
#TIMI

#This is the scripts to making each of the premeta files from the raw input files
rm(list=ls())
options(stringsAsFactors = FALSE)
library("data.table", lib.loc="/home/tj241/R_Packages")
source("/scratch/tj241/demo/Test/Manhattan_source.R")
source("/scratch/tj241/demo/Test/Manhattan/turboman_source.R")
####################
#UKBB - HAIL
####################
input=fread(paste0("/scratch/asb38/ukbiobank/CAD.INTERMEDIATE.allvariants.birthYear.sex.multiethnic.HRC.UK10K.1000G.no.varfilter.rel.excl.sel.3rd.deg.logreg.052318_chr1_22.tsv"), colClasses=c("character"))
input=as.data.frame(input)

temp=strsplit(input[,1], "[:]")
temp2=do.call(rbind, temp)

#Second allele is effect Allele
baseline=temp2[,3]
effect=temp2[,4]
CHR=as.numeric(temp2[,1])
POS=temp2[,2]

#MAF column is actually EAF
EAF=as.numeric(input[,5])
INFO=as.numeric(input[,4])

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(8,9,11)], EAF, INFO)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]
#There are perfect dupes with different RSIDs, fine to just drop one copy
output=unique(output)
#We have NA's, remove them
output=output[-which(is.na(output[,10])==T),]
#93095623 down to 93068690
#Filter INFO 0.4
output=output[as.numeric(output[,10])>0.4,]
#down to 68576130
#Filter MAF 0.00005
output=output[which(as.numeric(output[,9])>0.00005&as.numeric(output[,9])<0.99995),]
#down to 47902336
#Stitch on cases, 33206

Cases=rep(33206, length(output[,1]))
Effective_Cases=as.numeric(output[,10])*Cases
output=cbind(output, Cases, Effective_Cases)

#Resolve duplicates
dupes=output[output[,1]%in%output[duplicated(output[,1]),1],]
#These are all indels with flipped alleles which represent different variants. However, the only way to resolve them is to come up with a new naming system and apply it to every other study
#Decision: just drop all of these variants
output=output[-which(output[,1]%in%dupes[,1]),]
#down to 47900528

#We have 1044 with NA for betas/se/pval
output=output[-which(is.na(output[,6])==T|is.na(output[,7])==T|is.na(output[,8])==T),]

#47899484 left

#output[as.numeric(output[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output[,6]))
#summary(as.numeric(output[,7]))
#summary(as.numeric(output[,8]))

write.table(output, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/HAIL_adjusted_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")

####################
#UKBB - SAIGE
####################

options(scipen=999)
#Note, no scientific notation as they messed up their positions in the raw files
#Don't convert p-values, etc to numeric as they won't be in scientific notation

input=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/CAD_intermediate_GWAS_UK10K_allchr.SAIGE.txt", colClasses=c("character"))
input=as.data.frame(input)

#summary(pmin(nchar(input[,4]), nchar(input[,5])))
#temp=which(pmin(nchar(input[,4]), nchar(input[,5]))>1)
#Somehow we have complex indels
#ALLELE2 is the effect allele, the allele2af and allele2ac both refer to it so I can use as is
#Just leave complex indels in as is
baseline=toupper(input[,4])
effect=toupper(input[,5])
CHR=as.numeric(input[,1])
POS=as.numeric(input[,2])
#Position is numeric. However, no scientific notation so it's OK

#MAF column is actually EAF
EAF=as.numeric(input[,7])

#We need to get INFO from HAIL results, let's filter first


flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])
Cases=rep(33941, length(baseline))
N=rep(472335, length(baseline))

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(9,10,12)], EAF, Cases, N)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "Cases", "N")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]

#Must sort out duplicates before filters

print(dim(output))
output=unique(output)
print(dim(output))

#Now the flipped indels dupes
dupes=output[output[,1]%in%output[duplicated(output[,1]),1],]
print(dim(dupes))
output=output[-which(output[,1]%in%dupes[,1]),]
print(dim(output))

#Now I can apply filters
#Start with MAF

output=output[which(as.numeric(output[,9])>0.00005&as.numeric(output[,9])<0.99995),]
print(dim(output))
write.table(output, "/scratch/tj241/demo/Test/UKBB500K/CHD/SAIGE_500k_temporary.tsv", quote=F, row.names=F, col.names=T, sep="\t")
#dupes sorted and MAF filtered

#Time to append INFO
ref=fread(paste0("awk '{print $1,$4;}' /scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/CAD.INTERMEDIATE.allvariants.birthYear.sex.multiethnic.HRC.UK10K.1000G.no.varfilter.rel.excl.sel.3rd.deg.logreg.052318_chr1_22.tsv"), colClasses=c("character"))
ref=as.data.frame(ref)

temp=strsplit(ref[,1], "[:]")
temp2=do.call(rbind, temp)
baseline=temp2[,3]
effect=temp2[,4]
CHR=as.numeric(temp2[,1])
POS=as.numeric(temp2[,2])
#Again, pos is numeric since we set up no scientific notation earlier. 
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1
SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

ref[,1]=SNPID
colnames(ref)[1]="ALT_ID"
ref=unique(ref)


infopass=ref[which(as.numeric(ref[,2])>0.4),]
#infopass has duplicated IDs but these should be indels with flipped alleles, which have been removed entirely from output
output=output[output[,1]%in%infopass[,1],]
print(dim(output))

output2=merge(output, infopass, by="ALT_ID")
print(dim(output2))

finaloutput=output2[which(as.numeric(output2[,12])>0.4),]

Effective_Cases=as.numeric(finaloutput[,12])*as.numeric(finaloutput[,10])
finaloutput=cbind(finaloutput, Effective_Cases)

finaloutput=finaloutput[,c(1:9,12,10,13,11)]

print(head(finaloutput[as.numeric(finaloutput[,3])%%1000000==0,]))
#Check for NA beta/se/pval
print(summary(as.numeric(finaloutput[,6])))
print(summary(as.numeric(finaloutput[,7])))
print(summary(as.numeric(finaloutput[,8])))

write.table(finaloutput, "/scratch/tj241/demo/Test/UKBB500K/CHD/UKBB_SAIGE_premeta.tsv", quote=F, row.names=F, col.names=T, sep="\t")

#This takes an age to run, must submit it as job and leave it
#system(paste0("sbatch -N 1 -n 1 -c 10 -p long -J SAIGE -t 48:0:0 -o /scratch/tj241/demo/Test/Epic/UKBB_CHD_meta/Temp/saige.slurm /scratch/tj241/demo/Test/Epic/UKBB_CHD_meta/Temp/saige.sh"))


#####################
#Epic
####################
input=fread("/scratch/tj241/demo/Test/Epic/UKBB_CHD_meta/EPIC_CHD_gwas.tsv", colClasses=c("character"))
input=as.data.frame(input)
#input[as.numeric(input[,4])%%1000000==0,]
#this raw input is fine

files=c("/scratch/tj241/demo/Test/Epic/Incident_CHD/hce_cases", "/scratch/tj241/demo/Test/Epic/Incident_CHD/hce_noncases", "/scratch/tj241/demo/Test/Epic/Incident_CHD/quad_cases", "/scratch/tj241/demo/Test/Epic/Incident_CHD/quad_noncases")


temp1=fread(paste0(files[2], ".snpstats"))
temp1=as.data.frame(temp1)
temp2=fread(paste0(files[4], ".snpstats"))
temp2=as.data.frame(temp2)


discordance=cbind(temp1[,2], abs(as.numeric(temp1[,15])-as.numeric(temp2[,15])), temp1[,15], temp2[,15])
colnames(discordance)=c("rsid", "discordance","HCE_MAF", "QUAD_MAF")

output=merge(input, discordance, by="rsid", all.x=T, sort=F)
relative_discordance=as.numeric(output[,47])/as.numeric(output[,29])
output=cbind(output, relative_discordance)


ratio=pmin(as.numeric(output[,48])/as.numeric(output[,49]), as.numeric(output[,49])/as.numeric(output[,48]))
output=cbind(output, ratio)
#write.table(output, "/scratch/tj241/demo/Test/Epic/UKBB_CHD_meta/EPIC_CHD_gwas_full.tsv", quote=F, row.names=F, col.names=T, sep="\t")
#output=fread("/scratch/tj241/demo/Test/Epic/UKBB_CHD_meta/EPIC_CHD_gwas_full.tsv", colClasses=c("character"))
#output=as.data.frame(output)
#output[as.numeric(output[,4])%%1000000==0,]
#This is also correct
#In SNPtest, ALLELEB is the effect allele
baseline=output[,5]
effect=output[,6]
#no indels
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(output[,3], ":", output[,4], "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(output[which(flip==1),3], ":", output[which(flip==1),4], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

EAF=(2*as.numeric(output[,16])+as.numeric(output[,15]))/2/(as.numeric(output[,14])+as.numeric(output[,15])+as.numeric(output[,16]))

output2=cbind(SNPID, output[,3],output[,4], baseline, effect, output[,c(44,45,42,9,29,47,50,51)], EAF)
colnames(output2)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "INFO", "MAF", "Discordance", "Rel_Discordance", "Ratio", "EAF")

output3=output2[which(as.numeric(output2[,9])>0.7&as.numeric(output2[,10])>0.01&as.numeric(output2[,11])<0.05&as.numeric(output2[,12])<0.07&as.numeric(output2[,13])>0.9),]
#Final filters we actually applied
#INFO>0.7
#MAF>0.01
#Discordance<0.05
#Rel Discordance<0.07
#Ratio>0.9

#5.55M Variants remaining of the original 11.77M
#5703735 Variants actually, corrected now
#stitch on number of cases, 10365 for all SNPs
Cases=rep(10365, length(output3[,1]))
Effective_Cases=as.numeric(output3[,9])*Cases
N=rep(22430, length(output3[,1]))
output3=cbind(output3, Cases, Effective_Cases, N)
#output3[as.numeric(output3[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output3[,6]))
#summary(as.numeric(output3[,7]))
#summary(as.numeric(output3[,8]))
write.table(output3[,c(1:8,14, 9, 15, 16, 17)], "/scratch/tj241/demo/Test/UKBB500K/CHD/EPIC_premeta.tsv", quote=F, row.names=F, col.names=T, sep="\t")


####################
#Cardiogram - 1KG
####################

cardiogram=fread("/scratch/asb38/ukbiobank/results/cad_1000g_results_for_metal.txt", colClasses=c("character"))
cardiogram=as.data.frame(cardiogram)
cardiogram=cardiogram[order(as.numeric(cardiogram[,2]), as.numeric(cardiogram[,3])),]


allelemapping=fread("/scratch/asb38/ukbiobank/ALL_1000G_phase1integrated_v3.legend", colClasses=c("character"))
allelemapping=as.data.frame(allelemapping)
#Pretty sure first allele is "other allele"
colnames(allelemapping)=c("CHR", "SNPID", "POS", "other_allele", "reference_allele")

snps=c("A","T","G","C")
nonindels=cardiogram[-which(cardiogram[,4]=="D"|cardiogram[,5]=="D"),]
indels=cardiogram[which(cardiogram[,4]=="D"|cardiogram[,5]=="D"),]

baseline=nonindels[,5]
effect=nonindels[,4]
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1
CHR=nonindels[,2]
POS=nonindels[,3]
SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])
nonindels=cbind(SNPID, nonindels)

baseline=allelemapping[,4]
effect=allelemapping[,5]
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1
CHR=allelemapping[,1]
POS=allelemapping[,3]
SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])
allelemapping=cbind(SNPID, allelemapping)

temp=substr(allelemapping[,3],1,5)
fails=allelemapping[which(temp=="MERGE"),]

allelemapping3=allelemapping[-which(allelemapping[,3]%in%fails[,3]),]
colnames(allelemapping3)[3]="RSID"
#all cardiogram IDs are unique
#all allelemapping3 IDs are unique
finalnonindels=merge(nonindels, allelemapping3, by="SNPID")
#all matched and found

#We only need to fix the indels

#So we fix indels and just stick it at the end of nonindels then reorder
newid=paste0(indels[,2], ":", indels[,3])
indels=cbind(newid, indels)

newid=paste0(allelemapping[,2], ":", allelemapping[,4])
allelemapping=cbind(newid, allelemapping)
#trim allelemapping down to indels as well
allelemapping2=allelemapping[-which(allelemapping[,6]%in%snps&allelemapping[,7]%in%snps),]
#write.table(allelemapping2, "/scratch/tj241/demo/Test/UKBB500K/CHD/mapping_indels.tsv", quote=F, row.names=F, col.names=T, sep="\t")

#remove variants which begins with MERGED_DEL_...
#Plan is to remove all of them and then merge whatever I can
temp=substr(allelemapping2[,4],1,5)
fails=allelemapping2[which(temp=="MERGE"),]

allelemapping3=allelemapping2[-which(allelemapping2[,4]%in%fails[,4]),]

#Just merge it in

finalindel=merge(indels, allelemapping3, by="newid")
#lose 1534 indels
#Miraculously there are no "-" notation SNPs anymore!! Let's keep them this way for now then.
finalnonindels=finalnonindels[,2:18]
finalindel=finalindel[,c(2:13,15:19)]
colnames(finalindel)[14]="RSID"


#Process nonindels first
baseline=finalnonindels[,5]
effect=finalnonindels[,4]
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1
CHR=finalnonindels[,2]
POS=finalnonindels[,3]
SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])
finalnonindels=cbind(SNPID, finalnonindels)

#write.table(finalnonindels, "/scratch/tj241/demo/Test/UKBB500K/CHD/Cardiogram_nonindels.tsv", quote=F, row.names=F, col.names=T, sep="\t")

#Now to process the indels
#all of them have del as 1 char
baseline=finalindel[,5]
effect=finalindel[,4]

flip=nchar(finalindel[,16])==1
SNPID=paste0(finalindel[,2], ":", finalindel[,3], "_", finalindel[,17], "_", finalindel[,16])
SNPID[which(flip==T)]=paste0(finalindel[which(flip==T),2], ":", finalindel[which(flip==T),3], "_", finalindel[which(flip==T),16], "_", finalindel[which(flip==T),17])
finalindel=cbind(SNPID, finalindel)

deletion=pmin(finalindel[,17], finalindel[,18])
insertion=pmax(finalindel[,17], finalindel[,18])

newbaseline=deletion
newbaseline[which(finalindel[,6]=="I")]=insertion[which(finalindel[,6]=="I")]
neweffect=deletion
neweffect[which(finalindel[,5]=="I")]=insertion[which(finalindel[,5]=="I")]

finalindel[,6]=newbaseline
finalindel[,5]=neweffect

final=rbind(finalnonindels, finalindel)
final=final[order(final[,3], final[,4]),]
#write.table(final, "/scratch/tj241/demo/Test/UKBB500K/CHD/Cardiogram_full.tsv", quote=F, row.names=F, col.names=T, sep="\t")

#final=fread("/scratch/tj241/demo/Test/UKBB500K/CHD/Cardiogram_full.tsv", colClasses=c("character"))
#final=as.data.frame(final)
mapping=fread("/scratch/asb38/ukbiobank/snp.counts.txt1")
mapping=as.data.frame(mapping)

#apparantly all variants' RSIDs were found in the mapping file
colnames(mapping)[1]="RSID"

finaloutput=merge(final, mapping, by="RSID")

finaloutput2=finaloutput[,c(2,4,5,7,6,10:12,8,9,19)]
Effective_Cases=as.numeric(finaloutput[,9])*as.numeric(finaloutput[,19])
finaloutput2=cbind(finaloutput2, Effective_Cases)
colnames(finaloutput2)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO", "Cases", "Effective_Cases")
finaloutput2=finaloutput2[order(as.numeric(finaloutput2[,2]), as.numeric(finaloutput2[,3])),]

write.table(finaloutput2, "/scratch/tj241/demo/Test/UKBB500K/CHD/Cardiogram_premeta.tsv", quote=F, row.names=F, col.names=T, sep="\t")

####################
#Cardiogram - Exome
####################

input=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/MICAD.NA.EUR.ExA.Consortium.logistic.201404231.FINAL.txt", colClasses=c("character"))
input=as.data.frame(input)
missing=input[is.na(input[,23])==T|is.na(input[,22])==T,]
praveen1=fread("/scratch/backup/bioinformatics/illumina_manifests/exome_plus/version_2/5_other_output/version_1-0_restore/3_summary_files/v1_0_restore_info.csv", colClasses=c("character"))
praveen1=as.data.frame(praveen1)
colnames(praveen1)[4]="MarkerName"
test1=merge(missing, praveen1, by="MarkerName")
#perfect mapping file, but not for indel alleles
notmissing=input[-which(is.na(input[,23])==T|is.na(input[,22])==T),]
temp=test1[,c(1,26,27)]
tempoutput=merge(missing, temp,by="MarkerName")
tempoutput[,22]=tempoutput[,25]
tempoutput[,23]=tempoutput[,26]

#chr/pos sorted
output=rbind(notmissing, tempoutput[,1:24])
#write.table(output, "/scratch/tj241/demo/Test/UKBB500K/CHD/exome_temp_chrpos.tsv", quote=F, row.names=F, col.names=T, sep="\t")
indels=output[output[,2]=="d",]
#Pull chr/pos from 1kg files

nonindels=output[-which(output[,2]=="d"),]
test2=merge(indels, praveen1, by="MarkerName")

#Attempt to extract indel info from the strands
temp=gsub("(.*)(.)(\\[)(\\-)(\\/)(.+)(])(.*)", "\\2\\3\\4\\5\\6\\7", test2[,35])
#This is the correct one to use as I checked all 4 and they only differ for the last indel which was manually checked vs ukbb

test2[,2]=gsub("(.*)(.)(\\[)(\\-)(\\/)(.+)(])(.*)", "\\2", test2[,35])
test2[,3]=gsub("(.*)(.)(\\[)(\\-)(\\/)(.+)(])(.*)", "\\2\\6", test2[,35])
test2[,23]=as.numeric(test2[,23])-1
#Must make numeric!!
#Check perfect round positions at the end!!
exomefull=rbind(nonindels, test2[,1:24])
#This had better work
input=exomefull
#drop nonautosomal
input=input[-which(input[,22]%in%c("M", "X", "XY")),]
CHR=input[,22]
POS=input[,23]

baseline=toupper(input[,3])
effect=toupper(input[,2])

#Allele1 is assumed to be the effect allele until I get confirmation from Adam
#Confirmed
EAF=as.numeric(input[,4]) 
Cases=as.numeric(input[,20])
N=as.numeric(input[,21])+Cases
INFO=rep(1, length(Cases))
Effective_Cases=Cases

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(8:10)], EAF,INFO, Cases,Effective_Cases, N)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO","Cases","Effective_Cases", "N")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]
#There are perfect dupes with different RSIDs, fine to just drop one copy
output=unique(output)
#184249

#Hold up, there's duplicated variants here, same variant but different sample size/beta/se/pval/MAF. Need to filter down to unique set.
#we drop one with the fewest cases
dupes=output[duplicated(output[,1])==T, 1]

nondupes=output[-which(output[,1]%in%dupes),]
duplicates=output[which(output[,1]%in%dupes),]

dupesoutput=NULL
removal=NULL
for(i in 1:length(unique(duplicates[,1])))
{
    temp=duplicates[duplicates[,1]==unique(duplicates[,1])[i],]
    if(length(temp[,1])!=2)print(paste0("Different Number: ",length(temp[,1])))
    #all are pairs
    if(temp[1,10]==temp[2,10]&temp[1,11]==temp[2,11])
    {
        print(paste0("Same cases and controls: ",i))
        removal=c(removal, i)
    }
    temp=temp[order(-as.numeric(temp[,10]), -as.numeric(temp[,11])),]
    tempoutput=temp[1,]
    dupesoutput=rbind(dupesoutput, tempoutput)
}
#This still doesnt work, we have duplicated variants with the same chr/pos/alleles/cases/controls but somehow different MAF/effect/pval
#Just drop them
dupesoutput=dupesoutput[-removal,]
#29 lost

finaloutput=rbind(nondupes, dupesoutput)
finaloutput=finaloutput[order(as.numeric(finaloutput[,2]), as.numeric(finaloutput[,3])),]
#183644

#finaloutput[as.numeric(finaloutput[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(finaloutput[,6]))
#summary(as.numeric(finaloutput[,7]))
#summary(as.numeric(finaloutput[,8]))
write.table(finaloutput, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/Exome_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")


####################
#Cardiogram - Metabo
####################

input=fread("/scratch/curated_genetic_data/association_results/chd/EUR_ancestry_metabochip_79K_phenoscanner_april2017.tab", colClasses=c("character"))
input=as.data.frame(input)

CHR=input[,1]
POS=input[,3]

baseline=toupper(input[,6])
effect=toupper(input[,5])
#No indels
#Allele1 is effect allele
EAF=as.numeric(input[,7])
INFO=rep(1,length(effect))

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

totalsamples=as.numeric(input[,19])

#Getting cases too hard, we never got it for a couple of studies
#We will apply constant ratio of case:control given total sample size
#case:control = 73856:147712 from spreadsheet

Cases=totalsamples*73856/(147712+73856)
Effective_Cases=Cases

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(11:13)], EAF, INFO, Cases, Effective_Cases, totalsamples)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO","Cases","Effective_Cases", "N")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]
#There are perfect dupes with different RSIDs, fine to just drop one copy
output=unique(output)
#79037

#output[as.numeric(output[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output[,6]))
#summary(as.numeric(output[,7]))
#summary(as.numeric(output[,8]))

write.table(output, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/Metabo_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")


####################
#GerMIFs1-7 - not 3
####################

input1=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GerMIFSI.EUR.CAD.06072017.LZ.tab", colClasses=c("character"))
input1=as.data.frame(input1)
input2=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GerMIFSII.EUR.CAD.06072017.LZ.tab", colClasses=c("character"))
input2=as.data.frame(input2) 
input3=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GerMIFSIV.EUR.CAD.06072017.LZ.tab", colClasses=c("character"))
input3=as.data.frame(input3)
input4=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GerMIFSV.EUR.CAD.06072017.LZ.tab", colClasses=c("character"))
input4=as.data.frame(input4)
input5=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GerMIFSVI.EUR.CAD.12122017.LZ.tab", colClasses=c("character"))
input5=as.data.frame(input5)
input6=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GerMIFSVII.EUR.CAD.13062018.LL.tab", colClasses=c("character"))
input6=as.data.frame(input6)

#GerMIFS6 has problems, need an incl list
#incllist=read.table("/scratch/tj241/demo/Test/UKBB500K/CHD/GerMIFS6_incl_list.txt", header=F)
#keeplist is 2830433
#only found 2796168 of the 5108393 total 
#input2=input2[which(input2[,1]%in%incllist[,1]),]
#Leave GerMIFs6 as is and just flag all gwsig vars into a watch list, also add the chr1 false positive from EPIC to that watch list

summary(as.factor(input1[,6]))
summary(as.factor(input1[,7]))
summary(as.factor(input2[,6]))
summary(as.factor(input2[,7]))
summary(as.factor(input3[,6]))
summary(as.factor(input3[,7]))
summary(as.factor(input4[,6]))
summary(as.factor(input4[,7]))
summary(as.factor(input5[,6]))
summary(as.factor(input5[,7]))
summary(as.factor(input6[,6]))
summary(as.factor(input6[,7]))
#No indels on any of them
input=NULL
input[[1]]=input1
input[[2]]=input2
input[[3]]=NULL
input[[4]]=input3
input[[5]]=input4
input[[6]]=input5
input[[7]]=input6

numbers=c(622,1188,"NA",940,2392, 1639, 3062)
totals=c(2143,2426,"NA",2068,3929, 2825, 6524)

for(i in c(1,2,4:7))
{
    baseline=input[[i]][,6]
    effect=input[[i]][,7]
    #no indels
    flip=rep(0, length(baseline))
    flip[which(effect<baseline)]=1

    SNPID=paste0(input[[i]][,4], ":", input[[i]][,5], "_", baseline, "_", effect)
    SNPID[which(flip==1)]=paste0(input[[i]][which(flip==1),4], ":", input[[i]][which(flip==1),5], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

    #B is the effect allele
    EAF=(2*as.numeric(input[[i]][,10])+as.numeric(input[[i]][,9]))/2/(as.numeric(input[[i]][,8])+as.numeric(input[[i]][,9])+as.numeric(input[[i]][,10]))
    output=cbind(SNPID, input[[i]][,4],input[[i]][,5], baseline, effect, input[[i]][,c(14:17)], EAF)
    colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "INFO", "EAF")

    output2=output[which(output[,9]>0.4),]
                #INFO filter 0.4
                Cases=rep(numbers[i], length(output2[,1]))
                Effective_Cases=as.numeric(output2[,9])*as.numeric(Cases)
                N=rep(totals[i], length(output2[,1]))
                output2=cbind(output2, Cases, Effective_Cases, N)
                
                output2=output2[order(as.numeric(output2[,2]), as.numeric(output2[,3])),]
                
                #output2[as.numeric(output2[,3])%%1000000==0,]
                #Check for NA beta/se/pval
                #summary(as.numeric(output2[,6]))
                #summary(as.numeric(output2[,7]))
                #summary(as.numeric(output2[,8]))
                
                #summary(as.numeric(output2[,10]))
                #summary(as.numeric(output2[,9]))
                
    write.table(output2[,c(1:8,10,9,11,12,13)], paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/GerMIFS",i,"_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")

}


#####################
#deCODE
#####################
#Now for deCODE
#There are complex indels, cannot do anything
#Just process everything for now unless decode gets back to us
#WARNING! COMPLEX INDELS ARE CODED ALPHABETICALLY! Consequently they are not necessarily shorter indel first. 
input=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/DECODE.Coronary_Artery_Disease_All_06112017.18122017.txt", colClasses=c("character"))
input=as.data.frame(input)

input[,4]=substr(input[,4],4,nchar(input[,4]))
#Get rid of x chromosome
input=input[-which(input[,4]=="X"),]

baseline=input[,6]
effect=input[,7]
#Allele B is effect allele
EAF=as.numeric(input[,8])

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(input[,4], ":", input[,5], "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(input[which(flip==1),4], ":", input[which(flip==1),5], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, input[,4],input[,5], baseline, effect, input[,c(12:15)], EAF)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "INFO", "EAF")

finaloutput=output[which(output[,9]>0.4),]
#Stitch on cases, 38918
Cases=rep(38918, length(finaloutput[,1]))
Effective_Cases=as.numeric(finaloutput[,9])*Cases
N=rep(356875, length(finaloutput[,1]))
finaloutput=cbind(finaloutput, Cases, Effective_Cases, N)

#write.table(finaloutput, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/deCODE_unfiltered_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")

#Remove imperfect duplicates completely
finaloutput=finaloutput[-which(finaloutput[,1]%in%finaloutput[duplicated(finaloutput[,1]),1]),]
#3486 dropped

#They prefiltered by INFO score to 0.8
#Apply MAF filter 0.00005
#finaloutput=temp
finaloutput=finaloutput[which(as.numeric(finaloutput[,10])>0.00005&as.numeric(finaloutput[,10])<0.99995),]
#from 18398244 to 16595653

#finaloutput[as.numeric(finaloutput[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(finaloutput[,6]))
#summary(as.numeric(finaloutput[,7]))
#summary(as.numeric(finaloutput[,8]))

write.table(finaloutput[,c(1:8,10,9,11,12,13)], paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/deCODE_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")
#Note, complex indels not sorted but no other study have them so they're fine as is
#Note, 1k variants are duplicated perfectly but different EAF. This suggests different population but METAL automatically drops second one systematically. Must get confirmation from decode.
#Decision: just drop them since no response from deCODE


#####################
#Greek
#####################

#Now for Greek Dataset
input=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/GWAS_HRC.GREEK.DATASET.CAD.EUR.12122017.OG.txt", colClasses=c("character"))
input=as.data.frame(input)

#MAF=pmin((2*input[,6]+input[,7])/2/(input[,6]+input[,7]+input[,8]),(2*input[,8]+input[,7])/2/(input[,6]+input[,7]+input[,8]))
#Needs filtering. I have MAF calculated and it comes with INFO score
#Just INFO>0.4 no MAF filter
#No indels

baseline=input[,4]
effect=input[,5]
#B is effect allele
EAF=(2*as.numeric(input[,8])+as.numeric(input[,7]))/2/(as.numeric(input[,6])+as.numeric(input[,7])+as.numeric(input[,8]))

#no indels
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(input[,2], ":", input[,3], "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(input[which(flip==1),2], ":", input[which(flip==1),3], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, input[,2],input[,3], baseline, effect, input[,c(12:15)], EAF)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "INFO", "EAF")

output2=output[which(output[,9]>0.4&as.numeric(output[,10])<=0.95&as.numeric(output[,10])>=0.05),]
#Apply 0.4 INFO filter and 0.05 MAF filter
#stitch on number of cases, 335 for all SNPs
Cases=rep(335, length(output2[,1]))
Effective_Cases=as.numeric(output2[,9])*Cases
N=rep(1485, length(output2[,1]))
output2=cbind(output2, Cases, Effective_Cases, N)

#output2[as.numeric(output2[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output2[,6]))
#summary(as.numeric(output2[,7]))
#summary(as.numeric(output2[,8]))

write.table(output2[,c(1:8,10,9,11,12)], paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/Greek_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")


#########################
#HUNT - SAIGE
#########################
#We're happy with HUNT as well
input=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/CAD_HUNT.allchr.SAIGE.sorted.MAC.R2.MAF.hwe.dat", colClasses=c("character"))
input=as.data.frame(input)

colnames(input)=c("SNPID","CHR","POS","Allele0","Allele1","AC","AF","N","BETA","SE","Tstat","p.value","p.value.NA","Is.SPA.converge","varT","varTstar","MAC","R2","MAF","HWE_PVAL")

INFO=as.numeric(input[,18])
#Prefiltered to 0.3 INFO
#14614 have no INFO, presumably because they were genotyped like last time. Set their INFO to 1
INFO[which(is.na(INFO)==T)]=1

CHR=as.numeric(input[,2])
POS=input[,3]

baseline=toupper(input[,4])
effect=toupper(input[,5])
#ALLELE1 is taken to be the effect allele confirmed. AC and AF are all for Allele1.

#There are indels
#All indels are good, single variant deletion and no "-" so can process together
EAF=as.numeric(input[,7])

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(9,10,12)], EAF, INFO)
#output=output[order(as.numeric(output[,2]), as.numeric(output[,3]), output[,1]),]

colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO")
#Stitch on cases, 7710

Cases=rep(7710, length(output[,1]))
Effective_Cases=as.numeric(output[,10])*Cases
N=input[,8]
output2=cbind(output, Cases, Effective_Cases, N)

#INFO is prefiltered to 0.3
#We filter MAF to 0.0005

output2=output2[which(as.numeric(output2[,9])<=0.9995&as.numeric(output2[,9])>=0.0005),]
#Lose 8M variants, 16143079 left

#output2[as.numeric(output2[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output2[,6]))
#summary(as.numeric(output2[,7]))
#summary(as.numeric(output2[,8]))

output2=output2[order(as.numeric(output2[,2]), as.numeric(output2[,3]), output2[,1]),]

write.table(output2, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/HUNT_SAIGE_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")


#########################
#PARTNERS
#########################
#Process Partners
input=fread("/scratch/asb38/ukbiobank/CAD_meta_analysis_study_results/PHB_CAD_INTERMEDIATE_relevant_chips_ancestries_infoMAFcases_filtered_N_cases_N_total_meta_logORs_with_AVGINFO_0102191.tbl", colClasses=c("character"))
input=as.data.frame(input)

temp=strsplit(input[,1], "[:]")
temp2=do.call(rbind, temp)

CHR=temp2[,1]
POS=temp2[,2]

baseline=toupper(input[,3])
effect=toupper(input[,2])
#No indels
#Allele1 is effect allele
EAF=as.numeric(input[,4])
INFO=as.numeric(input[,20])
Cases=as.numeric(input[,17])
Effective_Cases=INFO*Cases
N=as.numeric(input[,18])

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(8:10)], EAF, INFO, Cases, Effective_Cases, N)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO", "Cases", "Effective_Cases", "N")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]
#output[as.numeric(output[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output[,6]))
#summary(as.numeric(output[,7]))
#summary(as.numeric(output[,8]))

write.table(output, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/PARTNERS_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")


################
#TIMI
################
#Process TIMI
input=fread("/scratch/asb38/TIMI/TIMI.EUR.CAD.age_sex_hxaf_hxt2d_batch_pc1_5.chr1_22.CR.FILT.20181102.txt", colClasses=c("character"))
input=as.data.frame(input)

#filter out non autosomals (bad chrs) and ones where chrs dont match
input=input[which(input[,22]%in%c(1:22)&input[,22]==input[,2]),]
#lost 35556
CHR=input[,22]
POS=input[,23]

baseline=toupper(input[,4])
effect=toupper(input[,6])

#We have indels but they have front allele. I can just proceed as expected.
#summary(pmin(nchar(baseline), nchar(effect)))


EAF=as.numeric(input[,25])
EAF[which(effect!=input[,26])]=1-as.numeric(input[which(effect!=input[,26]),25])
INFO=as.numeric(input[,17])
Cases=as.numeric(input[,27])
Effective_Cases=INFO*Cases
N=as.numeric(input[,8])

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(24, 10, 12)], EAF, INFO, Cases, Effective_Cases, N)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO", "Cases", "Effective_Cases", "N")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]

#drop all duplicate variants
output=output[-which(output[,1]%in%output[duplicated(output[,1]),1]),]
#3204 dupes, mostly indels but some are different positions in build 38 but same in hg19, e.g.chr1:121605507:C:T  and chr1:121853089:C:T 

#output[as.numeric(output[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output[,6]))
#summary(as.numeric(output[,7]))
#summary(as.numeric(output[,8]))


write.table(output, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/TIMI_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")

#########################
#HUNT - OLD
#########################
#We're happy with HUNT as well
input=fread("/scratch/asb38/ukbiobank/HUNT/CAD.intermediate.all.imputed.results_Cook.txt", colClasses=c("character"))
input=as.data.frame(input)
input=input[-which(input[,2]==23),]
#drop chr23, no other dataset has it

#info mapping file
info=fread("/scratch/asb38/ukbiobank/HUNT/CAD.intermediate.all.imputed.results.info.bed", colClasses=c("character"))
info=as.data.frame(info)

#same number of rows after removal of chr23 from input
INFO=as.numeric(info[,15])

CHR=info[,2]
POS=info[,3]
baseline=toupper(info[,5])
effect=toupper(info[,6])
flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1
SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

temp=cbind(SNPID, INFO)
#Problem, many variants don't have INFO score. Many are very rare but some are quite common
#Those are genotyped hence have info score of 1
temp[which(is.na(INFO)==T),2]=1
temp[,2]=as.numeric(temp[,2])

CHR=input[,2]
POS=input[,3]

baseline=toupper(input[,5])
effect=toupper(input[,6])
#ALLELE0 is taken to be the effect allele contrary to BOLT output manual. Adam gave the go ahead.
#This has been confirmed by the people who made the dataset
#There are indels
#All indels are good, single variant deletion and no "-" so can process together
EAF=as.numeric(input[,7])
#Actually EAF is just A1_freq despite the fact we want A0 Freq as A0 is our effect allele. Turns out, the freq column was also mislabelled.


flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(16,17,14)], EAF)
#output=output[order(as.numeric(output[,2]), as.numeric(output[,3]), output[,1]),]
output2=merge(output, temp, by="SNPID")
colnames(output2)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO")
#Stitch on cases, 7710

Cases=rep(7710, length(output2[,1]))
Effective_Cases=as.numeric(output2[,10])*Cases
output2=cbind(output2, Cases, Effective_Cases)

#INFO is prefiltered to 0.3
#We filter MAF to 0.0005

output2=output2[which(as.numeric(output2[,9])<=0.9995&as.numeric(output2[,9])>=0.0005),]
#Lose 8M variants, 16143117 left

#output2[as.numeric(output2[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output2[,6]))
#summary(as.numeric(output2[,7]))
#summary(as.numeric(output2[,8]))

output2=output2[order(as.numeric(output2[,2]), as.numeric(output2[,3]), output2[,1]),]

write.table(output2, paste0("/scratch/tj241/demo/Test/UKBB500K/CHD/HUNT_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")


#########################
#Biobank Japan
#########################
#Process Biobank Japan
input=fread("zcat /rds/project/jmmh2/rds-jmmh2-projects/coronary_genetics/million_hearts_project/Temp/bbj/CAD.auto.rsq07.mac10.txt.gz", colClasses=c("character"))
input=as.data.frame(input)


CHR=input[,1]
POS=input[,2]

baseline=toupper(input[,4])
effect=toupper(input[,5])
#Indels present but coded correctly, e.g. A/AAT
#Allele1 is effect allele
EAF=as.numeric(input[,7])
INFO=as.numeric(input[,19])
Cases=rep(29319 , length(input[,19]))
Effective_Cases=INFO*Cases
N=as.numeric(input[,8])

flip=rep(0, length(baseline))
flip[which(effect<baseline)]=1

SNPID=paste0(CHR, ":", POS, "_", baseline, "_", effect)
SNPID[which(flip==1)]=paste0(CHR[which(flip==1)], ":", POS[which(flip==1)], "_", effect[which(flip==1)], "_", baseline[which(flip==1)])

output=cbind(SNPID, CHR,POS, baseline, effect, input[,c(9,10,12)], EAF, INFO, Cases, Effective_Cases, N)
colnames(output)=c("ALT_ID", "CHR", "POS", "BASELINE", "EFFECT", "BETA", "SE", "PVAL", "EAF", "INFO", "Cases", "Effective_Cases", "N")
output=output[order(as.numeric(output[,2]), as.numeric(output[,3])),]

#75 indels with duplicated IDs. Removed them all
output2=output[-which(output[,1]%in%output[which(duplicated(output[,1])),1]),]
#output[as.numeric(output[,3])%%1000000==0,]
#Check for NA beta/se/pval
#summary(as.numeric(output[,6]))
#summary(as.numeric(output[,7]))
#summary(as.numeric(output[,8]))



write.table(output2, paste0("/rds/project/jmmh2/rds-jmmh2-projects/coronary_genetics/million_hearts_project/Temp/bbj/BBJ_premeta.tsv"), quote=F, row.names=F, col.names=T, sep="\t")

