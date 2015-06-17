#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("data.table")
library(phyloseq)
library(plyr)
library(ggplot2)
library(data.table)
setwd("/Users/adina/virome-paper/analysis-files")

==============
#LOAD FUNCTIONAL DATA
abund <- read.delim(sep='\t', file="./abundance_table.norm.txt",header=TRUE, strip.white=TRUE, row.names=1)

#Use the following to normalize the unnormalized abundance file
#d2 <- abund/rep(colSums(abund),each=length(abund$X77C_4.2_IND))
#abund <- d2
#write.table(x=d2, file="./abundance_table.norm.txt", sep="\t", col.names=TRUE, row.names=TRUE, quote=FALSE)
meta <- read.delim(sep=',', file="./metadata3.csv", header=TRUE, strip.white=TRUE, row.names=1)
meta <- data.frame(meta, cage_str = as.character(meta$cage))
ann <- read.delim(sep='\t', file="./giant-ann-fam.txt",header=FALSE, strip.white=TRUE, row.names=1)
colnames(ann) <- c("md5", "ident", "length", "V4", "V5", "V6","V7", "V8", "V9", "evalue", "bitscore", "ss", "ss_select", "l1", "l2", "l3", "t1", "t2","t3","t4","t5","t6", "t7")
ann_data_matrix <- as.matrix(ann)
annotation <- tax_table(ann_data_matrix)
abundance_data_norm_matrix <- as.matrix(abund)
abundance <- otu_table(abundance_data_norm_matrix, taxa_are_rows=TRUE)
metadata <- sample_data(meta)

#phyloseq object with annotations
all <- phyloseq(metadata, abundance, annotation)

#exploring temperate presence in virus fractions
all_temp_ind <- subset_samples(all, ext_frac == "ind")
all_temp_vlp <- subset_samples(all, ext_frac == "vlp")
mdf_ind <- psmelt(all_temp_ind)
mdf_ind_gt0 <- subset(mdf_ind, mdf_ind$Abundance > 0)
mdf_ind_gt0 <- subset(mdf_ind, mdf_ind$Abundance > 7.178e-5) #mean value of Abudances
mdf_vlp <- psmelt(all_temp_vlp)
mdf_vlp_gt0 <- subset(mdf_vlp, mdf_vlp$Abundance > 0)
mdf_vlp_gt0 <- subset(mdf_vlp, mdf_vlp$Abundance > 2.430e-6) #mean value of Abudances
t<-merge(mdf_ind_gt0, mdf_vlp_gt0, by="OTU")
t2<-unique(t$OTU)
blah <-ann[row.names(ann) %in% t2,]
summary(blah)
 
#exploring fatty acids
all_fatty <- subset_taxa(all, l3 == "Fatty Acids, Lipids, and Isoprenoids")
all_fatty <- subset_samples(all_fatty, ext_frac == "total")
all_fatty <- subset_samples(all_fatty, diet == "MF" | diet == "LF")
mdf <- psmelt(all_fatty)
mdf <- subset(mdf, l2 != "NA")
f <- ddply(mdf, .(cage, l3, ext_frac, diet, Sample, Age, l2), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(cage, l2, ext_frac, diet, Age), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
p = ggplot(f2, aes_string(x="Age", y="MEAN"))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p+xlab('')+theme_bw() +geom_point(stat="identity", size=3)+geom_errorbar(limits, width=0)+facet_grid(l2 ~ cage)# + theme(panel.grid.major=element_blank(), panel.grid.minor element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=10))
write.table(f, file="bacilli.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

summary(aov.result <- aov(SUM ~ l2*cage, data = f))
print(model.tables(aov.result,"means")) 
boxplot(SUM~l2*cage,data=f) 

TukeyHSD(aov.result)

l = length(levels(f$l1))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f, l1 == levels(f$l1)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			stat_sig[i] = as.character(f_dat$l1)
			print(TukeyHSD(aov.result))}
				}





#writing out subsets
vlp <- subset_samples(all, ext_frac == "ind" & cage_str == "50")
towrite <- otu_table(vlp)
write.table(towrite, file="si-50ind.txt", quote=FALSE, sep="\t", row.names=TRUE, col.names=TRUE)

#Exploring specific Taxa and functions represented
all2 <- subset_taxa(all, t3 == "Bacilli")
mdf <- psmelt(all2)
mdf <- subset(mdf, l3 == "Phages, Prophages, Transposable elements, Plasmids")
#mdf <- subset(mdf, l3 == "Carbohydrates")
#mdf <- subset(mdf, l3 == "Motility and Chemotaxis")
#summarizing average over experiment
f <- ddply(mdf, .(cage_str, t3, l3, ext_frac, diet2, Sample, Age), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(cage_str, t3, l3, ext_frac, diet2, Age), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
#gets rid of day 1 numbers, only want significance after 3 weeks on same diet
f3 <- subset(f2, diet2 != "MF-LF")
f3 <- subset(f3, diet2 != "LF-MF")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
p = ggplot(f3, aes_string(x="Age", y="MEAN"))
p+xlab('')+theme_bw() + facet_grid(ext_frac~cage_str)+theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", size=3)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=10))#+scale_y_log10()
write.table(f, file="bacilli.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)


#Overall characterizations
mdf <- psmelt(all)
f <- ddply(mdf, .(t2, cage_str, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(t2, cage_str, Age, ext_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- subset(f2, t2 != "NA")
f2 <- subset(f2, t2 != "")
f3 <- subset(f2, ext_frac == "vlp")
p = ggplot(f3, aes_string(x="t2", y="MEAN", color="t2"))
p+facet_grid(Age~cage_str)+xlab('')+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", size=5)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), legend.position="none")+ylab("Relative Abundance")+theme(text=element_text(size=20))+theme(axis.text.x = element_text(angle = 90, hjust=0, size=10))
#Standardized Bar charts (Figure 5 and 6)
p = ggplot(f3, aes_string(x="Age", y="MEAN", fill="t2"))
p+facet_grid(~cage_str)+xlab('')+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_bar(position="fill", stat="identity")+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), legend.position="none")+ylab("Relative Normalized Abundance")+theme(text=element_text(size=20))+theme(axis.text.x = element_text(angle = 90, hjust=0, size=10))

f <- ddply(mdf, .(l3, cage_str, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(l3, cage_str, Age, ext_frac), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- subset(f2, l3 != "NA")
f2 <- subset(f2, l3 != "")
f3 <- subset(f2, ext_frac == "vlp")
p = ggplot(f3, aes_string(x="l3", y="MEAN", color="l3"))
p+facet_grid(Age~cage_str)+xlab('')+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", size=5)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), legend.position="none")+ylab("Relative Abundance")+theme(text=element_text(size=20))+theme(axis.text.x = element_text(angle = 90, hjust=0, size=10))
#Standardized Bar charts (Figure 5 and 6)
p = ggplot(f3, aes_string(x="Age", y="MEAN", fill="l3"))
p+facet_grid(~cage_str)+xlab('')+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_bar(position="fill", stat="identity")+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank(), legend.position="none")+ylab("Relative Normalized Abundance")+theme(text=element_text(size=20))+theme(axis.text.x = element_text(angle = 90, hjust=0, size=10))

#to print out into tabular format for other plotting, we used PRISM for publication 
t <- subset(f2, ext_frac == "vlp")
t2<-t[order(-t$MEAN),]
t2$MEAN2 <- t2$MEAN*100
t2$SE2 <- t2$SE*100
write.table(t2, file="test.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#Specific subsystems 
mdf <- psmelt(all)
f <- ddply(mdf, .(l3, diet2, Age, ext_frac, Sample, cage_str), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(l3, diet2, Age, ext_frac, cage_str), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
f2 <- subset(f2, diet2 != "LF-MF")
f2 <- subset(f2, diet2 != "MF-LF")
subsys = "Motility and Chemotaxis"
f3 <- subset(f2, l3 == subsys)
f3 <- subset(f3, ext_frac == "vlp")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
temp = subset(f2, ext_frac == "total")
f3$diet2 = factor(f3$diet2, levels=c("LF1","MF*","LF2", "MF1","LF*","MF2"))
f2$l3 <- factor(f2$l3, levels=temp$l3[order(-temp$MEAN)])
p = ggplot(f3, aes_string(x="diet2", y="MEAN", shape="ext_frac", color="cage_str"))
p+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", size=7)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=25, family="Helvetica"), legend.position="none")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=15))+scale_colour_manual(values=c("orange","skyblue"))+xlab("")#+ggtitle(subsys)+xlab("")


#NMDS All fractions
all <- phyloseq(metadata, abundance)
#all <- phyloseq(metadata, abundance, annotation)
GPdist = phyloseq::distance(all, "bray")
GPNMDS = ordinate(all, "NMDS", GPdist)
#3/6 3/12 , 3/27, 4/2
#all_t <- subset_samples(all, date == "3/27/13" | date == "4/2/13")
#p2 = plot_ordination(all_t, GPNMDS, color="ext_frac", shape="cage_str")
p2 = plot_ordination(all, GPNMDS, color="ext_frac", shape="cage_str")+theme_bw()+theme(text=element_text(size=20))+geom_point(size=4)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p2 + guides(colour=guide_legend(title="DNA Fraction"), shape=guide_legend(title="Baseline diet"))
ggsave("FigureS6.eps")
dev.off()


p2 = plot_ordination(all, GPNMDS, shape="cage_str")
p2 + geom_point(size = 6)+ theme_bw() + theme(text = element_text(size=15))+theme_bw()+scale_x_continuous(limits = c(-.5, .75))+scale_y_continuous(limits = c(-.5, .5))+theme_bw()
x <- data.frame(GPNMDS$points)
y <- merge(x, meta, by="row.names")
p2 + geom_point(size = 5)+ theme_bw() + theme(text = element_text(size=15))+xlim(-.6000, .6000)+ylim(-.6000, .6000)+theme_bw()+ guides(colour=guide_legend(title="DNA Fraction"), shape=guide_legend(title="Baseline diet"))+scale_shape_manual("Baseline diet", values=c(16,17), labels=c("LF","MF"))
x <- data.frame(GPNMDS$points)
y <- merge(x, meta, by="row.names")
y <- data.frame(y, Age2 = as.character(y$Age))
meta <- data.frame(meta, Age = as.character(meta$Age))
y <- subset(y, ext_frac == "vlp")
ggsave("Fig7b.eps")
ggplot(data = y, aes(MDS1, MDS2)) +geom_point(data=subset(y, y$cage_str=="50"), colour="blue", size=5, aes(shape = Age2, size=10)) + geom_point(data=subset(y, y$cage_str=="77"), size=5, colour="orange", aes(shape=Age2, size=10))+scale_x_continuous(limits = c(-.3000,.6))+scale_y_continuous(limits = c(-.6,.4))+theme_bw()+theme(text=element_text(size=15)) +guides(shape=guide_legend(title="Day"))+scale_shape_manual("Day", values=c(16,1,17,15)) +theme(legend.key=element_rect(fill="white",colour="white"))+theme(legend.background = element_blank())
dev.off()
y <- subset(y, ext_frac == "vlp")
y <- subset(y, ext_frac == "ind")
ggsave("Fig7c.eps")
ggplot(data = y, aes(MDS1, MDS2)) +   geom_point(data=subset(y, y$cage_str=="50"), colour="blue", size=5, aes(shape = Age2, size=10)) + geom_point(data=subset(y, y$cage_str=="77"), size=5, colour="orange", aes(shape=Age2, size=10))+scale_x_continuous(limits = c(-.3000,.6))+scale_y_continuous(limits = c(-.6,.4))+theme_bw()+theme(text=element_text(size=15)) +guides(shape=guide_legend(title="Day"))+scale_shape_manual("Day", values=c(16,1,17,15))+theme(legend.key=element_rect(fill="white",colour="white"))+theme(legend.background = element_blank())
dev.off()
y <- subset(y, ext_frac == "total")
ggsave("Fig7a.eps")
ggplot(data = y, aes(MDS1, MDS2)) +   geom_point(data=subset(y, y$cage_str=="50"), colour="blue", size=5, aes(shape = Age2, size=10)) + geom_point(data=subset(y, y$cage_str=="77"), size=5, colour="orange", aes(shape=Age2, size=10))+scale_y_continuous(limits = c(-.1,.2))+scale_x_continuous(limits = c(-.6,-.3))+theme_bw()+theme(text=element_text(size=15)) +guides(shape=guide_legend(title="Day"))+scale_shape_manual("Day", values=c(16,17,15))+theme(legend.key=element_rect(fill="white",colour="white"))+theme(legend.background = element_blank())
dev.off()
adonis(GPdist ~ ext_frac, perm=9999, as(sample_data(all),"data.frame"))

#Generating NMDS Baseline distances
#x <- data.frame(GPNMDS$points)
#y <- merge(x, meta, by="row.names")
x <- read.delim(sep='\t', file="./nmds.txt",header=FALSE, strip.white=TRUE, row.names=1)
colnames(x)=c('dist','day')
#x$day <- as.character(x$day)
x2 <- merge(x, meta, by="row.names")
#x2 <- subset(x2, cage_str == "77")
f <- ddply(x2, .(day, ext_frac, cage_str), summarise, MEAN=mean(dist), SE=sd(dist)/sqrt(length(dist)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f$ext_frac <- factor(f$ext_frac, levels=c("total","ind","vlp"))
f <- subset(f, ext_frac == "total")
#f <- subset(f, day != "2")
f$day <- c("0","0","1","1","2","2")
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
c ="blue"
p = ggplot(f, aes_string(x="day",y="MEAN"),colour=c)
p + theme_bw()+geom_line(stat="identity", colour=c,aes(group=cage_str,linetype=cage_str), size=2)+geom_errorbar(limits, width=.1,colour=c)+geom_point(stat="identity",colour=c)+theme(axis.text.x = element_text(angle=90, hjust=1, size=5))+ylab("Average Distance from Baseline")+theme(text =element_text(size=15))

#Stats on NMDS
x3 <- subset(x2, cage_str == "77")
x2 <- x3
vlp = subset(x2, ext_frac == "ind")
vlp$day_c <- as.character(vlp$day)
summary(aov(dist~day_c, data=vlp))
aov.result<-aov(dist~day_c, data=vlp)
TukeyHSD(aov.result)
ind = subset(x2, ext_frac == "vlp")
ind$day_c <- as.character(ind$day)
summary(aov(dist~day_c, data=ind))
aov.result<-aov(dist~day_c, data=ind)
TukeyHSD(aov.result)



#Looking at clustering
veganotu = function(physeq) {
    require("vegan")
    OTU = otu_table(physeq)
    if (taxa_are_rows(OTU)) {
        OTU = t(OTU)
    }
    return(as(OTU, "matrix"))
}

#Examples on clustering
all <- phyloseq(metadata, abundance, annotation)
all <- phyloseq(metadata, abundance)
all_vlp <- all
all_vlp <- subset_samples(all, diet != 'LF-MF')
all_vlp <- subset_samples(all_vlp, diet != 'MF-LF')
all_vlp <- subset_samples(all, Age == "0" | Age == "43")
#all_vlp <- subset_samples(all_vlp, cage_str == "77")
#all_vlp <- subset_samples(all_vlp, diet2 != 'LF*')
#all_vlp <- subset_samples(all_vlp, diet2 != 'MF*')
meta <- sample_data(all_vlp)
v_all <- veganotu(all_vlp)
bd <- vegdist(v_all, method="bray")
adonis(bd ~ Age, perm=9999, as(sample_data(all_vlp),"data.frame"))
adonis(bd ~ ext_frac, perm=9999, as(sample_data(all_vlp),"data.frame"))
groups <- meta$Age
mod <- betadisper(bd, groups)
anova(mod)
permutest(mod, pairwise = TRUE, permutations=9999)
(mod.HSD <- TukeyHSD(mod))
clust <- hclust(bd, method="average")
setEPS()
postscript("suppfig4.eps")
plot(clust)
dev.off()

#Functions that are statistically significant
all <- phyloseq(metadata, abundance, annotation)
all_ss <- subset_samples(all, cage_str == "77")
all_ss <- subset_samples(all_ss, diet2 != "MF-LF")
all_ss <- subset_samples(all_ss, diet2 != "LF-MF")
mdf <- psmelt(all_ss)
f <- ddply(mdf, .(l3, diet2, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f3 <- subset(f, ext_frac == "ind") #change to vlp also
l = length(levels(f3$l3))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f3, l3 == levels(f3$l3)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			stat_sig[i] = as.character(f_dat$l3[1])
			print(TukeyHSD(aov.result))}
				}
				
#Contigs that are significant				
#mdf_phage = subset(mdf, l3 == "Phages, Prophages, Transposable elements, Plasmids")
#mdf_phage = subset(mdf, l3 == "Nitrogen Metabolism")
mdf_phage = subset(mdf, l3 == "Motility and Chemotaxis")
f <- ddply(mdf_phage, .(OTU, diet2, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f3 <- subset(f, ext_frac == "ind")
l = length(unique(f3$OTU))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f3, OTU == unique(f3$OTU)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			pvalue = summary(aov.result)[[1]][["Pr(>F)"]][1]
			stat_sig[i] = as.character(f_dat$OTU[1])
			print(TukeyHSD(aov.result))
			x<-c(unique(f_dat[1]),as.character(f_dat$SUM[1]), as.character(f_dat$SUM[2]),as.character(f_dat$SUM[3]),as.character(f_dat$SUM[4]), as.character(f_dat$SUM[5]),as.character(f_dat$SUM[6]), as.character(f_dat$SUM[7]),as.character(f_dat$SUM[8]), as.character(f_dat$SUM[9]), pvalue)
			write.table(x, file="ind.tsv", append=TRUE, quote=FALSE, sep="\t ", eol="\n", row.names = FALSE, col.names=FALSE)}
	}

#mdf_phage = subset(mdf, l3 == "Phages, Prophages, Transposable elements, Plasmids")
#mdf_phage = subset(mdf, l3 == "Nitrogen Metabolism")
mdf_phage = subset(mdf, l3 == "Motility and Chemotaxis")
f <- ddply(mdf_phage, .(OTU, diet2, Age, ext_frac, Sample), summarise, SUM=sum(Abundance))
f3 <- subset(f, ext_frac == "vlp")
l = length(unique(f3$OTU))
stat_sig = rep(0, l)
for (i in 1:l){
	f_dat <- subset(f3, OTU == unique(f3$OTU)[i])
	#print(f_dat)
	summary(aov.result <- aov(SUM ~ diet2, data = f_dat))
	if(summary(aov.result)[[1]][["Pr(>F)"]][1] <= 0.05 & summary(aov.result)[[1]][["Pr(>F)"]][1] != "NaN"){
			print(f_dat)
			print(summary(aov.result <- aov(SUM ~ diet2, data = f_dat)))
			pvalue = summary(aov.result)[[1]][["Pr(>F)"]][1]
			stat_sig[i] = as.character(f_dat$OTU[1])
			print(TukeyHSD(aov.result))
			x<-c(unique(f_dat[1]),as.character(f_dat$SUM[1]), as.character(f_dat$SUM[2]),as.character(f_dat$SUM[3]),as.character(f_dat$SUM[4]), as.character(f_dat$SUM[5]),as.character(f_dat$SUM[6]), as.character(f_dat$SUM[7]),as.character(f_dat$SUM[8]), as.character(f_dat$SUM[9]), pvalue)
			write.table(x, file="vlp.tsv", append=TRUE, quote=FALSE, sep="\t ", eol="\n", row.names = FALSE, col.names=FALSE)}
	}

	
# Exploring core proteins (Figure 8)
all <- phyloseq(metadata, abundance)
all_virus <- subset_samples(all, ext_frac != "total")
all_virus <- subset_samples(all, ext_frac != "total" & cage_str == "50")
abund <- otu_table(all_virus)
virus_core <- abund[apply(abund, MARGIN=1, function(x) all(x > 0)),]
#virus_core_phy <- phyloseq(virus_core, metadata, annotation)
virus_core_phy <- phyloseq(virus_core, metadata)
mdf <- psmelt(virus_core_phy)
f <- ddply(mdf, .(OTU, cage_str, ext_frac, Age), summarise, MEAN=mean(Abundance), SE=sd(Abundance)/sqrt(length(Abundance)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
f2 = subset(f, ext_frac == "ind")
#labels annotated
subsetx <- f2[which(f2$OTU == "10397_10396_IND" | f2$OTU == "10398_10397_IND" | f2$OTU == "11291_11290_IND" | f2$OTU == "11296_11295_IND" | f2$OTU == "11327_11326_IND" | f2$OTU == "11343_11342_IND" | f2$OTU == "1392833_97998_IND" | f2$OTU == "1541070_99692_IND" | f2$OTU == "1541404_99769_IND" | f2$OTU == "1858_1857_VLP" | f2$OTU == "2033_2032_VLP" | f2$OTU == "2250_2249_VLP" | f2$OTU == "22767_6492_TOT" | f2$OTU == "2366_2365_VLP" | f2$OTU == "2368_2367_VLP" | f2$OTU == "26_6575_TOT" | f2$OTU == "2650_1950_TOT" | f2$OTU == "2701_2700_VLP" | f2$OTU == "276638_14923_VLP" | f2$OTU == "28_3490_VLP" | f2$OTU == "3327_3326_VLP" | f2$OTU == "3349_3348_VLP" | f2$OTU == "4585_5883_TOT" | f2$OTU == "6916_6915_IND" | f2$OTU == "74376_10776_VLP" | f2$OTU == "9907_9906_IND"),]
temp <- subset(f2, cage_str == "50" & Age == "0")
temp_foo <- subset(subsetx, cage_str == "50" & Age == "0")
f2$OTU <- factor(f2$OTU, levels=temp$OTU[order(-temp$MEAN)])
subsetx$OTU <- factor(subsetx$OTU, levels=temp_foo$OTU[order(-temp_foo$MEAN)])
p = ggplot(f2, aes_string(x="OTU", y="MEAN"))
p+theme_bw()+geom_point(stackt="identity", aes(color=ext_frac))+geom_errorbar(limits, width=0)+facet_grid(Age~cage_str)

#clustering significant contigs
abund_sig <- abund[which(rownames(abund) ==  "11070_5415_VLP" |  rownames(abund) ==  "11392_11391_IND" |  rownames(abund) ==  "11495_11494_IND" |  rownames(abund) ==  "11500_11499_IND" |  rownames(abund) ==  "11573_11572_IND" |  rownames(abund) ==  "11577_11576_IND" |  rownames(abund) ==  "1168485_95035_IND" |  rownames(abund) ==  "11895_11894_IND" |  rownames(abund) ==  "1301_1300_VLP" |  rownames(abund) ==  "1395_1394_VLP" |  rownames(abund) ==  "1460_1459_VLP" |  rownames(abund) ==  "1541036_99683_IND" |  rownames(abund) ==  "1545310_100359_IND" |  rownames(abund) ==  "1546282_100483_IND" |  rownames(abund) ==  "1568425_102182_IND" |  rownames(abund) ==  "189_188_IND" |  rownames(abund) ==  "2052_2051_IND" |  rownames(abund) ==  "2281_2280_VLP" |  rownames(abund) ==  "268318_14671_VLP" |  rownames(abund) ==  "278888_15308_VLP" |  rownames(abund) ==  "279518_15392_VLP" |  rownames(abund) ==  "282627_15745_VLP" |  rownames(abund) ==  "3314_3313_VLP" |  rownames(abund) ==  "3398_3397_VLP" |  rownames(abund) ==  "345_344_IND" |  rownames(abund) ==  "47_46_IND" |  rownames(abund) ==  "485381_86166_IND" |  rownames(abund) ==  "6233_6232_IND" |  rownames(abund) ==  "81561_11054_VLP" |  rownames(abund) ==  "8182_8181_IND" |  rownames(abund) ==  "8233_8232_IND" |  rownames(abund) ==  "9169_9168_IND" |  rownames(abund) ==  "9619_9618_IND"),]
bcd <- vegdist(abund_sig, method="bray")
#bcd<-vegdist(foo, method="bray")
clust <- hclust(bcd, method="average")
plot(clust)
#saved cluster names to cluster.txt file
cluster <-read.delim(sep='\t', file="./cluster.txt",header=FALSE, strip.white=TRUE, row.names=1)
colnames(cluster) <- c("g")
cluster <- data.frame(cluster)
cluster$g <- as.character(cluster$g)
bcd <- vegdist(abund_sig, method="bray")
adonis(decostand(abund_sig,"total") ~ g,data=cluster,perm=9999)
groups <- cluster$g
mod <- betadisper(bcd, groups)
anova(mod)
permutest(mod, pairwise = TRUE, permutations=9999)
(mod.HSD <- TukeyHSD(mod))

#putting signif contigs into phyloseq object
abundance_data_norm_matrix <- as.matrix(abund_sig)
abundance <- otu_table(abundance_data_norm_matrix, taxa_are_rows=TRUE)
signif <- phyloseq(metadata, abundance, annotation)
signif_77 <- subset_samples(signif, cage_str == "77")
t<-merge_samples(signif_77, "Age")
d <- psmelt(t)
d$log_abund <- log2(d$Abundance)
ggplot(d, aes(x=Sample, y=OTU, fill=Abundance)) + geom_tile() + scale_fill_gradient(low="#000033", high="#FF3300")+theme_bw()+xlab('Experimental Day')+ylab('Assembled Contig ID')



#looking into 35 contigs
head(annotation)
ann_sig <- ann[which(rownames(ann) ==  "11070_5415_VLP" |  rownames(ann) ==  "11392_11391_IND" |  rownames(ann) ==  "11495_11494_IND" |  rownames(ann) ==  "11500_11499_IND" |  rownames(ann) ==  "11573_11572_IND" |  rownames(ann) ==  "11577_11576_IND" |  rownames(ann) ==  "1168485_95035_IND" |  rownames(ann) ==  "11895_11894_IND" |  rownames(ann) ==  "1301_1300_VLP" |  rownames(ann) ==  "1395_1394_VLP" |  rownames(ann) ==  "1460_1459_VLP" |  rownames(ann) ==  "1541036_99683_IND" |  rownames(ann) ==  "1545310_100359_IND" |  rownames(ann) ==  "1546282_100483_IND" |  rownames(ann) ==  "1568425_102182_IND" |  rownames(ann) ==  "189_188_IND" |  rownames(ann) ==  "2052_2051_IND" |  rownames(ann) ==  "2281_2280_VLP" |  rownames(ann) ==  "268318_14671_VLP" |  rownames(ann) ==  "278888_15308_VLP" |  rownames(ann) ==  "279518_15392_VLP" |  rownames(ann) ==  "282627_15745_VLP" |  rownames(ann) ==  "3314_3313_VLP" |  rownames(ann) ==  "3398_3397_VLP" |  rownames(ann) ==  "345_344_IND" |  rownames(ann) ==  "47_46_IND" |  rownames(ann) ==  "485381_86166_IND" |  rownames(ann) ==  "6233_6232_IND" |  rownames(ann) ==  "81561_11054_VLP" |  rownames(ann) ==  "8182_8181_IND" |  rownames(ann) ==  "8233_8232_IND" |  rownames(ann) ==  "9169_9168_IND" |  rownames(ann) ==  "9619_9618_IND"),]
summary(ann_sig)
write.table(ann_sig, file="annotation_of_sig.txt", quote=FALSE,sep="\t")
ann_sig_g1 <- ann[which(rownames(ann) ==  "11070_5415_VLP" |  rownames(ann) ==  "11495_11494_IND" |  rownames(ann) ==  "1168485_95035_IND" |  rownames(ann) ==  "11895_11894_IND" |  rownames(ann) ==  "1301_1300_VLP" |  rownames(ann) ==  "1541036_99683_IND" |  rownames(ann) ==  "1546282_100483_IND" |  rownames(ann) ==  "2052_2051_IND" |  rownames(ann) ==  "3314_3313_VLP" |  rownames(ann) ==  "345_344_IND" |  rownames(ann) ==  "47_46_IND" |  rownames(ann) ==  "485381_86166_IND" |  rownames(ann) ==  "8233_8232_IND" |  rownames(ann) ==  "9619_9618_IND"),]
summary(ann_sig_g1)
ann_sig_g2 <- ann[which(rownames(ann) ==  "11392_11391_IND" |  rownames(ann) ==  "11500_11499_IND" |  rownames(ann) ==  "11573_11572_IND" |  rownames(ann) ==  "11577_11576_IND" |  rownames(ann) ==  "1395_1394_VLP" |  rownames(ann) ==  "1460_1459_VLP" |  rownames(ann) ==  "1545310_100359_IND" |  rownames(ann) ==  "1568425_102182_IND" |  rownames(ann) ==  "189_188_IND" |  rownames(ann) ==  "2281_2280_VLP" |  rownames(ann) ==  "268318_14671_VLP" |  rownames(ann) ==  "278888_15308_VLP" |  rownames(ann) ==  "279518_15392_VLP" |  rownames(ann) ==  "282627_15745_VLP" |  rownames(ann) ==  "3398_3397_VLP" |  rownames(ann) ==  "6233_6232_IND" |  rownames(ann) ==  "81561_11054_VLP" |  rownames(ann) ==  "8182_8181_IND" |  rownames(ann) ==  "9169_9168_IND"),]
summary(ann_sig_g2)

ann_data_matrix_g1 <- as.matrix(ann_sig_g1)
annotation_g1 <- tax_table(ann_data_matrix_g1)
g1 <- phyloseq(annotation_g1, abundance, metadata)
mdf <- psmelt(g1)
f <- ddply(mdf, .(l3, diet2, Age, ext_frac, cage_str, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(l3, diet2, Age, ext_frac, cage_str), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
ggsave("group1.eps")
p = ggplot(f2, aes_string(x="Age", y="MEAN", shape="ext_frac", color="cage_str"))+facet_grid(ext_frac~cage_str)
p+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", sdize=7)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=25, family="Helvetica"), legend.position="none")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=15))+scale_colour_manual(values=c("blue","orange"))
dev.off()





ann_data_matrix_g2 <- as.matrix(ann_sig_g2)
annotation_g2 <- tax_table(ann_data_matrix_g2)
g2 <- phyloseq(annotation_g2, abundance, metadata)
mdf <- psmelt(g2)
f <- ddply(mdf, .(l3, diet2, Age, ext_frac, cage_str, Sample), summarise, SUM=sum(Abundance))
f2 <- ddply(f, .(l3, diet2, Age, ext_frac, cage_str), summarise, MEAN=mean(SUM), SE=sd(SUM)/sqrt(length(SUM)))
limits<-aes(ymin=MEAN-SE, ymax=MEAN+SE)
ggsave("group2.eps")
p = ggplot(f2, aes_string(x="Age", y="MEAN", shape="ext_frac", color="cage_str"))+facet_grid(ext_frac~cage_str)
p+theme_bw() + theme(strip.background = element_rect(colour='white',fill = 'white'))+geom_point(stat="identity", sdize=7)+geom_errorbar(limits, width=0)+ theme(panel.grid.major=element_blank(), panel.grid.minor = element_blank())+ylab("Relative Abundance")+theme(text=element_text(size=25, family="Helvetica"), legend.position="none")+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5,size=15))+scale_colour_manual(values=c("blue","orange"))
dev.off()
