library(corto)

############### creation of the NAT regulon #############################


#Normal adjacent tissue samples
#dataframes of samples
cts=read.table('input/RNA-seq_GNAME_VST_NAT.tsv', sep='\t', header = TRUE)
cts_nat<-cts[,-1]
rownames(cts_nat)<-cts[,1]

load("input/RNA-seq_GNAME_VST_NAT.rda")


#Tumor samples
cts_tumor=read.table('input/RNA-seq_GNAME_VST_TUM.tsv', sep='\t', header = TRUE)
cts_tum<-cts_tumor[,-1]
rownames(cts_tum)<-cts_tumor[,1]

load("input/RNA-seq_GNAME_VST_TUM.rda")


#centroids list
tfgenes=read.table('input/centroids_TFs.txt', header = TRUE)
tfgenes<-unlist(tfgenes)
tfgenes<-unname(tfgenes)

load("input/centroids_TFs.rda")
typeof(tfgenes)

### regulon


nthreads<-8
p<-1e-8
nbootstraps<-1000

load("output/regulon_TUMvsNAT/TUMvsNAT_reg_NAT-regulon_onlyDBCK.rda")

################### MRA analyisis #########################################

fname<-"output/mra_TUMvsNAT/mra_TUMvsNAT_reg_NAT-regulon_onlyDBCK.rda"

if(!file.exists(fname)){

    mra_TUMvsNAT_regNAT_regulon_dbck<-mra(cts_tum,cts_nat,regulon_dbck,nthreads=8,minsize=20)
    save(mra_TUMvsNAT_regNAT_regulon_dbck,file=fname)
} else {load(fname)}


############### plot of results ###########################################

png("plot/mra_TUMvsNAT_top10_onlyDBCK.png",w=600,h=700,res=100)
mraplot(mra_TUMvsNAT_regNAT_regulon_dbck,mrs=10)
dev.off()

png("plot/mra_TUMvsNAT_regNAT_regulon_dbck_todown.png",w=700,h=800,res=100)
toshow<-names(sort(mra_TUMvsNAT_regNAT_regulon_dbck$nes,decreasing=FALSE))
toshow

mraplot(mra_TUMvsNAT_regNAT_regulon_dbck,mrs=toshow)
dev.off()

png("plot/mra_TUMvsNAT_regNAT_regulon_dbck_topup.png",w=800,h=900,res=100)
toshow<-names(sort(mra_TUMvsNAT_regNAT_regulon_dbck$nes,decreasing=TRUE))
toshow

mraplot(mra_TUMvsNAT_regNAT_regulon_dbck,mrs=toshow)
dev.off()
