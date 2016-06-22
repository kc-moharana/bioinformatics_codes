
## SYSTEM REQUIREMENTS: UNIX, R , (OPTIONALLY R-Studio)
## OBJECTIVES: TO FIND Differentially Expressing GENEs  (DEGs) IN 4 PATHOLOGICAL CONDOTIONS IN <I>VITIS VINIFERA<\I>
##	BN INFECTED, VIRUS INFECTED, RECOVERD FROM INFECTION AND HEALTHY GRAPE LEAVES. 	



##DATA SOURCE: Vitis - transcription profiling by array 
#########################################################
#Total RNA was extracted from central leaf midribs and petioles from different V. vinifera cultivars in different conditions 
#(healthy, infected and recovered). Microarray analyses were conducted using different biological replicates for treatment.


# Nimblgen raw data should be in .pair file format; But NO RAW data was provided in GEO; So download series matrix file and extreact the expression data
# open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17894
# click : Series Matrix File(s) Link
# Or download fom here: ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE52nnn/GSE52540/matrix/GSE52540_series_matrix.txt.gz
# Extreact zipped file GSE52540_series_matrix.txt.gz file to GSE52540_series_matrix.txt
# On ubuntu terminal use command: 
# UNIX CMD: gunzip GSE52540_series_matrix.txt.gz


##CREATE EXPRESSION FILE FROM SERIES MATRIX FILE
##################################
#GSE52540_series_matrix.txt is a big file with all different information. Expression data will be present after the commant lines line ("!series_matrix_table_begin")
# Extract lines not starting with '!' symbols.
# using perl script to filter file:  perl -ne 'unless(/^!/) { print "$_"}' GSE52540_series_matrix.txt > GSE52540_series_matrix.filtered.txt
# number of probes found: 29973
# Number of columns: 31 (first col: probe id; second onwords: expressions), 30 samples;

##CREATE SAMPLE PHENOTYPE FILE FROM SERIES MATRIX FILE
###################################
# Extreact the comment lines from series file and open it in Excel file
# select line number 31 and 40, !Sample_geo_accession and !Sample_characteristics_ch1 (with healthy information and cultival info), respectively. Paste-special --> transpose.
#Separate Cultivar, Disease condition and replicates, by replacing commas with tabs
# Add Column headings: geo_ac	cultivar	status	replicate	condition
# See file named: sample_category.txt



##################################################################
##ANNOTATING PROBE SET IDS WITH V.vinifera GENE/TRANSCRIPT/PEPTIDE NAMES; and A. thaliana ORTHOLOGS
########################################################
# open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17894
# At the bottom, Download GPL17894_090918_Vitus_exp.ndf
# as we dont have any file describing gene locus information for each probe set id, we need to peform blastn with vitis transcriptome sequences.
# Extract Probe sequences from "GPL17894_090918_Vitus_exp.ndf"  file using Unix command cut to fetch 5th and 6th column, seq_id and probe_seq respectively. Later converting two columns in to two lines; first one seq_id and second sequence
# UNIX CMD:   cut -f 5,6  GPL17894_090918_Vitus_exp.ndf | perl -lane 'print ">$F[0]\n$F[1]" if $F[0] ne "EMPTY"' > GPL17894_090918_Vitus_exp.fa

# Download V. vinifera annoations from Phytozome. (https://phytozome.jgi.doe.gov). You may need to log-in as a user to download these annotations.
# log-in Go to speciese--> Select V. vinifera--> Bulk Data --> select "Vvinifera_145_Genoscope.12X.transcript_primaryTranscriptOnly.fa " and "Vvinifera_145_Genoscope.12X.annotation_info.txt"
# Click download selected files

# Using BLASTn to match Probe sequences to find corresponding Gene name. 
# Create blast database
#makeblastdb -in Vvinifera_145_Genoscope.12X.transcript_primaryTranscriptOnly.fa  -dbtype nucl
# Run Blastn
# UNIX CMD: blastn -max_target_seqs 1 -num_threads 3 -ungapped -outfmt 7 -word_size 22 -query GPL17894_090918_Vitus_exp.fa  -db Vvinifera_145_Genoscope.12X.transcript_primaryTranscriptOnly.fa -out  GPL17894_090918_Vitus_exp_probe_blast.txt

# Fetch the first Gene hit; Unicx command Grep is used
# UNIX CMD: grep -v "#" GPL17894_090918_Vitus_exp_probe_blast.txt | perl -lane 'if($h{$F[0]}){next; }$h{$F[0]}=1; print "$F[0]\t$F[1]"  ' > GPL17894_090918_Vitus_exp_probe_blast.tab


## Maping with V. vinifera gene annotation. using perl
# UNIX CMD: perl -lane 'if(scalar @F > 2 ){$h{$F[2]}=join("\t",@F); }else {print "$F[0]\t$F[1]\t$h{$F[1]}" if $h{$F[1]}; }  ' Vvinifera_145_Genoscope.12X.annotation_info.txt GPL17894_090918_Vitus_exp_probe_blast.tab > GPL17894_090918_Vitus_exp_probe_annotation_mapped.txt

## The annotation file contains A. tahliana gene anootations; useful for many useful annotation.
### create a small file containing only Probe_set_id and Gene names



###############################################################################################
## All data stored in "~/STUDY_material/BIOINFO/Microarray folder" 



###############################################################################################
# Alternative way // Automatics download and analysis using GEOquery; Needs active internet
# Upto this step can be done through GEOquery alsoe using the following command
# library(GEOquery)
# library(limma)
# gse_list = getGEO("GSE53179")
# gse = gse_list[["GSE53179_series_matrix.txt.gz"]]
# head(pData(gse))
###############################################################################################
 
 
##Installing required libraries
 source("http://bioconductor.org/biocLite.R")
 biocLite("affy")
 biocLite("limma") 

 
 
##Set working directory; thhis will be different for you
setwd("~/STUDY_material/BIOINFO/Microarray")

## reading expression data in to R object
exprs_dat <- read.table("GSE52540_series_matrix.filtered.txt", sep="\t", header = TRUE, check.names = FALSE, na.strings=c("NA", "-", "?") )
# little dataformatting
row.names(exprs_dat) <- exprs_dat$ID_REF
exprs_dat<- exprs_dat[,2:31]
head(exprs_dat)
#checking number of rows and columns
dim(exprs_dat)
#[1] 29971    30
#Convert to matrix
exprs_dat<- data.matrix(exprs_dat, rownames.force = NA)


##Reading phenotypic data/ disease condition into an R object.
samples <- read.table("sample_category.txt", sep="\t", header = TRUE, check.names = FALSE)
row.names(samples) = samples$geo_ac
#samples = samples[,2:5]


## NORMALIZATION
################################
## chks names are sorted or not
colnames(exprs_dat)
#[1] "GSM1269432" "GSM1269433" "GSM1269434" "GSM1269435" "GSM1269436" "GSM1269438" "GSM1269439" "GSM1269440" "GSM1269441" "GSM1269442" "GSM1269443" "GSM1269444"
#[13] "GSM1269445" "GSM1269446" "GSM1269447" "GSM1269448" "GSM1269449" "GSM1269450" "GSM1269451" "GSM1269452" "GSM1269453" "GSM1269454" "GSM1269455" "GSM1269456"
#[25] "GSM1269457" "GSM1269458" "GSM1269459" "GSM1269460" "GSM1269461" "GSM1269462"
samples$geo_ac
#[1] GSM1269432 GSM1269433 GSM1269434 GSM1269435 GSM1269436 GSM1269438 GSM1269439 GSM1269440 GSM1269441 GSM1269442 GSM1269443 GSM1269444 GSM1269445 GSM1269446 GSM1269447
#[16] GSM1269448 GSM1269449 GSM1269450 GSM1269451 GSM1269452 GSM1269453 GSM1269454 GSM1269455 GSM1269456 GSM1269457 GSM1269458 GSM1269459 GSM1269460 GSM1269461 GSM1269462
#30 Levels: GSM1269432 GSM1269433 GSM1269434 GSM1269435 GSM1269436 GSM1269438 GSM1269439 GSM1269440 GSM1269441 GSM1269442 GSM1269443 GSM1269444 GSM1269445 ... GSM1269462
## these data are already RMA normalized; 

# check distributiono of expression values
pdf("Box-plot.pdf")
 par(mar=c(7,2,1,2))
 boxplot(log2(exprs_dat),cex=0.5, las=2, col=sapply(samples$status,function(x) ifelse(x=="control","red",ifelse(x=="infected","orange","green"))),main="Distribution of expression values")
 legend("bottomright",inset=.02, legend=c("Healthy", "BN/virus infected", "Recovered"),  col=c("red", "orange", "green"), pch=15,  cex=0.8)
dev.off()
 
##Principal compnent analysisn using prcomp(); 
# if PCA is required for samples then, they must be the row names.
pc=prcomp(t(exprs_dat),center=T,scale.=T)	##using transformation t() function to create sample X gene table
pdf("PCA_plot.pdf")
par(mar=c(20,4,6,4))
par(mfrow=c(1,3)) 
plot(pc$x[,1:2],col=sapply(samples$status,function(x) ifelse(x=="control","red",ifelse(x=="infected","orange","green"))), pch=20)
# Some useful R tricks are commented 
#    text(pc$x[,1:2], labels =row.names(pc$x), cex=0.6, pos=2)
plot(pc$x[,2:3],col=sapply(samples$status,function(x) ifelse(x=="control","red",ifelse(x=="infected","orange","green"))), pch=20)
plot(pc$x[,1:3],col=sapply(samples$status,function(x) ifelse(x=="control","red",ifelse(x=="infected","orange","green"))), pch=20)
legend("topleft",inset=.02, legend=c("Healthy", "BN/virus infected", "Recovered"),  col=c("red", "orange", "green"), pch=20,  cex=0.5 )
dev.off()



##Limma package requires expression data to be in 'Expression set' objects (eset object)
# Creating ESET object from scart expression data
# SOURCE: http://stackoverflow.com/questions/25581769/creating-eset-object-from-preprocessed-expression-matrix

#load Biobase lib.
library("Biobase")

#annoationobject for exprs data using sample object
pd <- new("AnnotatedDataFrame", data = samples)


#Creating an ExpressionSet from pd and exprs_dat
eset <- ExpressionSet(assayData =exprs_dat, phenoData = pd )
#Check names
sampleNames(eset)[1:5]
# check conditions
eset$condition == "healthy"
#[1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[29] FALSE FALSE
eset$condition == "BN infected"
#[1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE FALSE FALSE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
#[29] FALSE FALSE

##  SIGNIFICANCE TEST USING LIMMA
#################################
library("limma")
# WE HAVE three main GROUPS NAMELY 'control' , 'infected' AND 'recovered'
# control=9, infected=17 AND recovered=4; Total =30
# design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3)))

## WE HAVE three five small GROUPS NAMELY "healthy", "BNinfected","VIRinfected","BN_VIRinfected", "recovered"
# 5 condition series; c(1,1,1,1,1,1,1,1,1,2,2,2,3,3,2,2,2,2,4,4,4,3,3,3,3,3,5,5,5,5)

#design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,1,1,1,2,2,2,3,3,2,2,2,2,4,4,4,3,3,3,3,3,5,5,5,5)))
#--
## Use which oone is useful


#ADD COLUMN NAMES
colnames(design) <- c("healthy", "BNinfected","VIRinfected","BN_VIRinfected", "recovered")
design

###############################################################################################
# COMPARING DATASETS
#-------------------------
#IN CONTRAST MATRIX WE HAVE TO SPECIFY WHAT PAIR WE WANT TO COMPARE, 

#===========================
#Analysis -1               #
#===========================
## How many genes change from healthy to BN infectin
## creating contrast matrix for BNinfected-healthy (change here to compare other combinations if you wise)
design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,1,1,1,2,2,2,3,3,2,2,2,2,4,4,4,3,3,3,3,3,5,5,5,5)))
#ADD COLUMN NAMES
colnames(design) <- c("healthy", "BNinfected","VIRinfected","BN_VIRinfected", "recovered")
cont.matrix <- makeContrasts(infVScontrol=BNinfected-healthy, levels=design)

##CREATE MArrayLM OBJECTS BY lmFit AND eBayes
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

##FIND AND FILTER THE OPGENES IN THE MATRIX ON DIFFERENT PARAMETERS 
# coef= column number or column name specifying which coefficient or #contrast of the linear model is of interest.
# adjust.method="none", "BH", "BY" and "holm"
# resort.by="logFC", "AveExpr", "t", "P", "p", "B" or "none"
# lfc=minimum absolute log2-fold-change required

TopGenes_all=topTable(fit2,number=Inf, adjust.method="BH",sort.by="B", p.value=1, lfc=0)
# write.csv(TopGenes_all, file="infected-control_allProbes_logFold.csv")

##EXPORT ALL PROBE SET IDS
#============================

# You may open this output file in excel and apply various filters on logFC and BH columns.
# Same thing can be done using R also: using various combination of pvalue and logFC cutoffs
# Use ONLY oNE suitable combination
# p.value <0.05 and log fold change (lfc) >= +-2
# TopGenes_DEGs=topTable(fit2, number=Inf, adjust.method="BH",sort.by="B", p.value=0.05, lfc=2)
# p.value <0.01 and log fold change (lfc) >= +-2
# TopGenes_DEGs=topTable(fit2, number=Inf, adjust.method="BH",sort.by="B", p.value=0.01, lfc=2)
# p.value <0.05
# TopGenes_DEGs=topTable(fit2, number=Inf, adjust.method="BH",sort.by="B", p.value=0.05, lfc=2)
# log fold change (lfc) >= +-2
# TopGenes_DEGs=topTable(fit2, number=Inf, adjust.method="BH",sort.by="B", p.value=1, lfc=2)

## Exporting top 2000 significant probe sets (adj.P.Val<0.05); applying BH correction
TopGenes_DEGs_BN_H=topTable(fit2, number=2000, adjust.method="BH",sort.by="logFC", p.value=0.05 )
write.csv(TopGenes_DEGs_BN_H, file="BNinfected-control_DEG_Probes_logFold.csv")

#VOLCANO PLOTS TO DISPALY DEGs
################################
pdf("Vitis_DEG_volcano_plot_With_genesBN_vs_healthy.pdf", width = 6, height=6)
par(mgp = c(1.5,0.4,0))
# TopGenes_all, plot all genes logFC vs adj.Pval in black dots
with(TopGenes_all, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot: Bois Noir disease vs Control (Pval<0.05 and logFC>2 ) ", xlab = "log2 fold change", ylab = "-log10(P-value)", ylim=c(0,14) ))
# TopGenes_DEGs, plot all genes logFC vs adj.Pval in red dots
with(TopGenes_DEGs_BN_H, points(logFC, -log10(adj.P.Val), pch=20, col = "red"))
# TopGenes_DEGs, label significant genes; be careful and change as per your need 
#with(subset(TopGenes_DEGs_BN_H, adj.P.Val <1e-10 & abs(logFC)>3), text(logFC,-log10(adj.P.Val), labels=row.names(TopGenes_DEGs_BN_H),  cex=.5,pos = 3 ))
legend("topright", inset=0.002, title = "Significantly DEGs (P.val<0.05)", legend=c("Yes", "no"), col = c("red", "black"), pch=20)
dev.off()

#HEATMAP PLOT
#######################
# sub set expression data for GEGs and plot heatmap
# read annotatoin_ file to print gene names on het map (See below how to cerate a thei file GPL17894_090918_Vitus_exp_probe_blast.tab)
grape_gene_annot= read.csv("GSE52540_RAW/GPL17894_090918_Vitus_exp_probe_blast.tab", sep="\t", header=F)
pdf("Vitis_DEG_HeatMap_With_genesBN_vs_healthy.pdf", width = 6, height=6)
row.names(grape_gene_annot)= grape_gene_annot$V1
heatmap(exprs_dat[row.names(topTable(fit2, number=100, adjust.method="BH",sort.by="logFC", p.value=0.05, lfc=2) ),], cexRow =0.5, labRow =grape_gene_annot[row.names(topTable(fit2, number=100, adjust.method="BH",sort.by="logFC", p.value=0.05, lfc=2) ),]$V2, main = "Top 100 DEGs in Bois Noir disease vs Control "  )
dev.off()



#===========================
#Analysis -2               #
#===========================
#
#  How many genes change from healthy to Virus infectin
# creating contrast matrix for VIRinfected-healthy (change here to compare other combinations if you wise)
cont.matrix <- makeContrasts(infVScontrol=VIRinfected-healthy, levels=design)

fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

TopGenes_all=topTable(fit2,number=Inf, adjust.method="BH",sort.by="logFC", p.value=1, lfc=0)
TopGenes_DEGs_VIR_H=topTable(fit2, number=2000, adjust.method="BH",sort.by="logFC", p.value=0.05 )
write.csv(TopGenes_DEGs_VIR_H, file="VIRinfected-control_DEG_Probes_logFold.csv")


## VOLCANO PLOTS TO DISPALY DEGs

pdf("Vitis_DEG_volcano_plot_With_genesVir_vs_healthy.pdf", width = 6, height=6)
par(mgp = c(1.5,0.4,0))

with(TopGenes_all, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot: Virus infected vs Control (Pval<0.05 and logFC>2 ) ", xlab = "log2 fold change", ylab = "-log10(P-value)", ylim=c(0,14) ))

with(TopGenes_DEGs_VIR_H, points(logFC, -log10(adj.P.Val), pch=20, col = "red"))
# TopGenes_DEGs, label significant genes; be careful and change as per your need 
legend("topright", inset=0.002, title = "Significantly DEGs (P.val<0.05)", legend=c("Yes", "no"), col = c("red", "black"), pch=20)
dev.off()
## HEATMAP PLOT

## sub set expression data for GEGs and plot heatmap
## read annotatoin_ file to print gene names on het map (See below how to cerate a thei file GPL17894_090918_Vitus_exp_probe_blast.tab)
#grape_gene_annot= read.csv("GSE52540_RAW/GPL17894_090918_Vitus_exp_probe_blast.tab", sep="\t", header=F)
pdf("Vitis_DEG_HeatMap_With_genesVIR_vs_healthy.pdf", width = 6, height=6)
row.names(grape_gene_annot)= grape_gene_annot$V1
heatmap(exprs_dat[row.names(topTable(fit2, number=100, adjust.method="BH",sort.by="logFC", p.value=0.05, lfc=2) ),], cexRow =0.5, labRow =grape_gene_annot[row.names(topTable(fit2, number=100, adjust.method="BH",sort.by="logFC", p.value=0.05, lfc=2) ),]$V2, main = "Top 100 DEGs in Virus infected vs Control "  )
dev.off()




#===========================
#Analysis -3               #
#===========================
## How many genes change from Infectio to recovery
## creating contrast matrix for recoverd-infectd  (change here to compare other combinations if you wise)
## create new factor ; merge all type of infection to infection ()
design <- model.matrix(~ 0+factor(c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3)))
colnames(design) <- c("healthy", "infected", "recovered")
cont.matrix <- makeContrasts(recovVSinfect=recovered-infected, levels=design)

fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

TopGenes_all=topTable(fit2,number=Inf, adjust.method="BH",sort.by="B", p.value=1, lfc=0)
TopGenes_DEGs_reco_inf=topTable(fit2, number=2000, adjust.method="BH",sort.by="logFC", p.value=0.05 )
write.csv(TopGenes_DEGs_reco_inf, file="Recoverd-infectio_DEG_Probes_logFold.csv")


#VOLCANO PLOTS TO DISPALY DEGs

pdf("Vitis_DEG_volcano_plot_With_genesrecovered_vs_infected.pdf", width = 6, height=6)
par(mgp = c(1.5,0.4,0))

with(TopGenes_all, plot(logFC, -log10(adj.P.Val), pch=20, main="Volcano plot: Recovered vs Phy+Infected (Pval<0.05 and logFC>2 ) ", xlab = "log2 fold change", ylab = "-log10(P-value)", ylim=c(0,14) ))

with(TopGenes_DEGs_reco_inf, points(logFC, -log10(adj.P.Val), pch=20, col = "red"))
# TopGenes_DEGs, label significant genes; be careful and change as per your need 
legend("topright", inset=0.002, title = "Significantly DEGs (P.val<0.05)", legend=c("Yes", "no"), col = c("red", "black"), pch=20)
dev.off()
#HEATMAP PLOT

## sub set expression data for GEGs and plot heatmap
## read annotatoin_ file to print gene names on het map (See above paragraphs for  how to cerate a the file GPL17894_090918_Vitus_exp_probe_blast.tab)
#grape_gene_annot= read.csv("GSE52540_RAW/GPL17894_090918_Vitus_exp_probe_blast.tab", sep="\t", header=F)
pdf("Vitis_DEG_HeatMap_With_genesRecoveredVsInfected.pdf", width = 6, height=6)
row.names(grape_gene_annot)= grape_gene_annot$V1
heatmap(exprs_dat[row.names(topTable(fit2, number=100, adjust.method="BH",sort.by="logFC", p.value=0.05, lfc=2) ),], cexRow =0.5, labRow =grape_gene_annot[row.names(topTable(fit2, number=100, adjust.method="BH",sort.by="logFC", p.value=0.05, lfc=2) ),]$V2, main = "Top 100 DEGs in Recovered vs Phy+infected "  )
dev.off()



## open all .csv files in excel and may use vlookup() function to map to desired Probe data annottion




