#OPLSDA for CON-LVD using strain-adjusted values based on y~ diet+strain 
#Note that although transformation was performed up until line 163, data were not transformed for strain-correction and OPLSDA

#setwd("C:\\Users\\UNC NRI\\Desktop\\Jing NRI\\Metabolomics\\G0_081418_8strains3diets_JX\\OPLSDA_usingStrAdj")

which(bad.LVD)
#These are 4 metabolites with NA or identical values across all samples in the CON-LVD dataset
#m.4-vinylphenol sulfate              m.daidzein             m.genistein             m.glycitein 
#165                     268                     318                     339

library(car)
library(lmtest)

residual<-matrix(data=NA,nrow=n.LVD, ncol=n.met)
colnames(residual)=temp.name
bp.pval<-rep(NA,length=n.met)
sw.pval<-rep(NA,length=n.met)
pb <- winProgressBar(title="R Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i]) next
        residual[,i]=resid(lm(meta.LVD[,i]~diet.LVD + strain.LVD))
        bp.pval[i]=bptest(lm(meta.LVD[,i]~diet.LVD + strain.LVD))$p.value
        sw.pval[i]=shapiro.test(residual[,i])$p.value
        
        info<-sprintf("%d%% done", round(i/n.met*100))
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum(bp.pval>=0.05,na.rm=TRUE) #482 homostasciticity 
sum(bp.pval<0.05,na.rm=TRUE) #166 heterostasciticity 
sum(sw.pval>=0.05,na.rm=TRUE) #160 normal 
sum(sw.pval<0.05,na.rm=TRUE) #488 not normally distributed 

sum((bp.pval>=0.05)&(sw.pval>=0.05),na.rm=TRUE) #111 meet both requirements without transformation

######(3-d) BoxCox transformation########
library(MASS)
meta.boxcox=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.boxcox)=temp.name
#http://rcompanion.org/handbook/I_12.html# 
pb <- winProgressBar(title="Transformation Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i])next
        model=lm(meta.LVD[,i]~diet.LVD + strain.LVD)
        Box=boxcox(model,lambda = seq(-6,6,0.1))
        Cox=data.frame(Box$x,Box$y) 
        Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
        lambda= Cox2[1, "Box.x"]# Extract the lambda with greatest log-likelihood 
        if(lambda==0){meta.boxcox[,i]=log(meta.LVD[,i])}else{
                meta.boxcox[,i]=(meta.LVD[,i] ^ lambda - 1)/lambda
        }
        
        rm(Box,Cox,Cox2,lambda)# or rm(list=c('Box','Cox'))
        
        info<-sprintf("%d th metabolite transformed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

res.boxcox<-matrix(data=NA,nrow= n.LVD, ncol=n.met)
colnames(res.boxcox)=temp.name
bp.boxcox<-rep(NA,length=n.met)
sw.boxcox<-rep(NA,length=n.met)

pb <- winProgressBar(title="Regression Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i])next
        res.boxcox[,i]=resid(lm(meta.boxcox[,i]~diet.LVD + strain.LVD))
        bp.boxcox[i]=bptest(lm(meta.boxcox[,i]~diet.LVD + strain.LVD))$p.value
        sw.boxcox[i]=shapiro.test(res.boxcox[,i])$p.value
        
        info<-sprintf("%d th metabolite done", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum((bp.boxcox>=0.05)&(sw.boxcox>=0.05),na.rm=TRUE) #449/652 meet both requirements after BoxCox
sum(((bp.boxcox>=0.05)&(sw.boxcox>=0.05))|((bp.pval>=0.05)&(sw.pval>=0.05))
    ,na.rm=TRUE) #461/652 meet both requirements before or after BoxCox

######(3-c) Tukey transformation########
meta.tukey=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.tukey)=temp.name
res.tukey<-matrix(data=NA,nrow= n.LVD, ncol=n.met)
colnames(res.tukey)=temp.name
bp.tukey<-rep(NA,length=n.met)
sw.tukey<-rep(NA,length=n.met)

library(rcompanion)#http://rcompanion.org/handbook/I_12.html 
pb <- winProgressBar(title="Transformation Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i])next
        meta.tukey[,i]=transformTukey(meta.LVD[,i],plotit=FALSE,verbose = FALSE)
        
        info<-sprintf("%d th metabolite transformed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

for(i in 1:n.met){
        if(bad.LVD[i]) next
        res.tukey[,i]=resid(lm(meta.tukey[,i]~diet.LVD + strain.LVD))
        bp.tukey[i]=bptest(lm(meta.tukey[,i]~diet.LVD + strain.LVD))$p.value
        sw.tukey[i]=shapiro.test(res.tukey[,i])$p.value
}
sum((bp.tukey>=0.05)&(sw.tukey>=0.05),na.rm=TRUE) #417/652 meet both requirements after Tukey transformation
sum(((bp.boxcox>=0.05)&(sw.boxcox>=0.05))|((bp.pval>=0.05)&(sw.pval>=0.05))
    |((bp.tukey>=0.05)&(sw.tukey>=0.05))
    ,na.rm=TRUE) #475/652 meet both requirements before or after BoxCox or Tukey transformation

######(3-d) Log transformation########
meta.logtr=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.logtr)=temp.name
res.logtr<-matrix(data=NA,nrow= n.LVD, ncol=n.met)
colnames(res.logtr)=temp.name
bp.logtr<-rep(NA,length=n.met)
sw.logtr<-rep(NA,length=n.met)

for (i in 1:n.met){
        if(bad.LVD[i])next
        meta.logtr[,i]=log(meta.LVD[,i])
        res.logtr[,i]=resid(lm(meta.logtr[,i]~diet.LVD + strain.LVD))
        bp.logtr[i]=bptest(lm(meta.log.LVD[,i]~diet.LVD + strain.LVD))$p.value
        sw.logtr[i]=shapiro.test(res.logtr[,i])$p.value
}
sum((bp.logtr>=0.05)&(sw.logtr>=0.05),na.rm=TRUE) #315/652 meet both requirements after Tukey transformation
sum(((bp.boxcox>=0.05)&(sw.boxcox>=0.05))|((bp.pval>=0.05)&(sw.pval>=0.05))
    |((bp.tukey>=0.05)&(sw.tukey>=0.05))|((bp.logtr>=0.05)&(sw.logtr>=0.05))
    ,na.rm=TRUE) #483/652 meet both requirements before or after BoxCox or Tukey or log transformation

######(3-e) combine transformation files########
Pval_assump<-rbind(bp.pval,sw.pval,bp.boxcox,sw.boxcox,bp.tukey,sw.tukey,bp.logtr,sw.logtr)
colnames(Pval_assump)=temp.name
rawGood<-apply(Pval_assump[1:2,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
boxGood<-apply(Pval_assump[3:4,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
tukGood<-apply(Pval_assump[5:6,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
logGood<-apply(Pval_assump[7:8,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))

apply(Pval_assump,1,function(x) sum(is.na(x)))
#bp.pval   sw.pval bp.boxcox sw.boxcox  bp.tukey  sw.tukey  bp.logtr  sw.logtr 
#4         4       4         4          4         4         4         4 

Good=rawGood|logGood|boxGood|tukGood
sum(Good,na.rm = TRUE) #483
sum(!Good,na.rm = TRUE) #165
sum(is.na(Good)) #4
Compliance<-rbind(Good,rawGood, boxGood, tukGood,logGood)
colnames(Compliance)=temp.name

#priority high-low: raw, boxcox, tukey, log
data.choice=vector(length=n.met)
for (i in 1:n.met){
        if(bad.LVD[i]){data.choice[i]="noChoice"}
        else if(rawGood[i]==TRUE) {data.choice[i]="raw"}
        else if(boxGood[i]==TRUE) {data.choice[i]="boxcox"}
        else if(tukGood[i]==TRUE) {data.choice[i]="tuk"}
        else if(logGood[i]==TRUE) {data.choice[i]="log"}
        else {data.choice[i]="noChoice"}
}
table(data.choice)
#boxcox      log noChoice      raw      tuk 
#350        8      169      111       14 
table(logGood[data.choice=="tuk"])#3/14 Tukey transformed metabolites can use log transformation to meet both assumptions
bp=barplot(table(data.choice),main="transformation for OLS-linear regression", 
           xlab="type of transformation", ylab="number of metabolites")
text(bp,15,table(data.choice), cex=1, pos=3)






######(3-f) correct metabolite values by strain########
#Data not transformed in anyway########################
#CON and line 1a were set as baseline (value=0) in regression
m=lm(meta.LVD[,1]~diet.LVD+strain.LVD)
names(m$coefficients)
#[1] "(Intercept)"           "diet.LVDLow Vitamin D" "strain.LVD1b"          "strain.LVD2a"          "strain.LVD2b"         
#[6] "strain.LVD3a"          "strain.LVD3b"          "strain.LVD9a"          "strain.LVD9b"

#strain adjusted response = y - beta2*X_strain
get_str.adjY=function(Y, diet, strain){
        model=lm(Y~diet+strain)
        Y.adj=Y
        len=length(Y)
        for (i in 1:len){
                str.i=as.character(strain[i])
                if (str.i!="1a"){
                        coef=model$coefficients
                        n=names(coef)
                        beta=coef[which(substr(n,nchar(n)-1,nchar(n))==str.i)] 
                        Y.adj[i]=Y[i]-beta  
                }
        }
        Y.adj
}
#test run
get_str.adjY(meta.LVD[,1], diet.LVD,strain.LVD)

meta.LVD.adj=meta.LVD
pb <- winProgressBar(title="Adjust for strain", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i]) next
        meta.LVD.adj[,i]=get_str.adjY(meta.LVD[,i], diet.LVD, strain.LVD)
        
        info<-sprintf("%d%% done", round(i/n.met*100))
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

library(plyr)
#Use plyr revalue() to change "Low Vitamin D" to "LVD"
d.LVD=revalue(diet.LVD, c("Low Vitamin D"="LVD"))
diet_str_LVD=paste(d.LVD, strain.LVD, sep = "-")
row.names(meta.LVD.adj)=diet_str_LVD
row.names(meta.LVD)=diet_str_LVD


#####OPLSDA################################################################################################################
library("ropls")
#http://huboqiang.cn/2016/03/03/RscatterPlotPCA#35-save-in-pdf-file   Helpful link for plotting
library(ggplot2)
library(grid)
library(gridExtra)

##replace NA values to 1 in order to run ropls package
rmNA=function(metafile){
        n=nrow(metafile)
        new=metafile
        uniq=apply(metafile,2,function(x) length(unique(x)))
        x=which(uniq==1)
        m=length(x)
        replace=rep(1,n)
        for (i in 1:m){
                new[,x[i]]=replace
        }
        new
}
meta.new.LVD=rmNA(meta.LVD)
meta.new.LVD.adj=rmNA(meta.LVD.adj)

pca.LVD.adj<-opls(meta.new.LVD.adj)
plot(pca.LVD.adj, typeVc = "x-score",parAsColFcVn = diet.LVD, parEllipsesL = TRUE)


oplsda_LVD=opls(meta.new.LVD, diet.LVD,predI = 1, orthoI = 1, permI=1000)#1000 permutations
#OPLS-DA
#72 samples x 648 variables and 1 response
#standard scaling of predictors and response(s)
#4 excluded variables (near zero variance)
#R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort  pR2Y   pQ2
#Total    0.131    0.607 -0.0999  0.32   1   1 0.061 0.111

oplsda_LVD.adj=opls(meta.new.LVD.adj, diet.LVD,predI = 1, orthoI = 1, permI=1000)#1000 permutations
#OPLS-DA
#72 samples x 648 variables and 1 response
#standard scaling of predictors and response(s)
#4 excluded variables (near zero variance)
#R2X(cum) R2Y(cum) Q2(cum) RMSEE pre ort  pR2Y   pQ2
#Total     0.11    0.748   0.266 0.256   1   1 0.001 0.001

oVIP.LVD.adj<-getVipVn(oplsda_LVD.adj, orthoL = TRUE) #648 metabolites, removed the 4.
pVIP.LVD.adj<-getVipVn(oplsda_LVD.adj, orthoL = FALSE) #648 metabolites, removed the 4.

LVD.oplscore.padj<-getScoreMN(oplsda_LVD.adj, orthoL=FALSE)
LVD.oplscore.oadj<-getScoreMN(oplsda_LVD.adj, orthoL=TRUE)
LVD.oplscore.adj=cbind(LVD.oplscore.padj,LVD.oplscore.oadj)
LVD.oplscore.adj=as.data.frame(LVD.oplscore.adj)
opl.perc.LVD.adj=c("p1 (4%)", "to1")
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA)
             ,panel.grid.major = element_blank(),panel.grid.minor = element_blank()
             ,strip.background=element_blank(),axis.text.x=element_text(colour="black")
             ,axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black")
             ,plot.margin=unit(c(1,1,1,1),"line"))
p<-ggplot(LVD.oplscore.adj,aes(x=p1,y=o1,color=diet.LVD))
p<-p+geom_point(size=3, alpha=0.5)+theme+xlab(opl.perc.LVD.adj[1]) + ylab(opl.perc.LVD.adj[2]) +
        geom_text(aes(label=factor(substr(unlist(data.LVD$`CC line_newname` ),4,5))),
                  hjust=-0.2, vjust=-0.2,size=3,fontface='bold') +
        ggtitle("Q2=0.266, pQ2=0.001")
p
ggsave("OPLSDA_LVD_strAdjust_dot.pdf", p,width=2, height=1.6, units="in", scale=3)

##add VIP score to metabolite list
addVIP=function(meta, VIP){
        n=ncol(meta)
        metalist=colnames(meta)
        VIP.new=rep(0, n)
        names(VIP.new)=metalist
        for (i in 1:n){
                current=metalist[i]
                position=which(names(VIP)==current)
                if(length(position)>0){
                        VIP.new[i]=VIP[position]
                }
                
                
        }
        VIP.new
}
orthVIP.new.LVD.adj=addVIP(meta.new.LVD.adj,oVIP.LVD.adj)
predVIP.new.LVD.adj=addVIP(meta.new.LVD.adj,pVIP.LVD.adj)

sum(orthVIP.new.LVD.adj>1.5, na.rm = TRUE) #85 metabolites with adj.orthVIP>1.5
sum(predVIP.new.LVD.adj>1.5, na.rm = TRUE) #86 metabolites with adj.predVIP>1.5

rmsVIP=sqrt(orthVIP.new.LVD.adj^2+predVIP.new.LVD.adj^2)/sqrt(2)

VIPvector=cbind(orthVIP.new.LVD.adj,predVIP.new.LVD.adj,rmsVIP)
write.table(VIPvector,"Output//LVD_strAdj_VIPmatrix.txt",sep = "\t",quote = F,col.names = NA)

library("ggplot2")
# Calculating Pearson's product-moment correlation
pdf("Output//VIPcorrelation.pdf")
rmsVIPcolor=ifelse(rmsVIP>1.5,'>1.5',ifelse(rmsVIP<1,'<1','1-1.5'))

p=ggplot(data.frame(VIPvector),aes(orthVIP.new.LVD.adj,predVIP.new.LVD.adj))
p+geom_point(aes(orthVIP.new.LVD.adj,predVIP.new.LVD.adj,color=rmsVIPcolor))+
        stat_smooth(method="lm",formula=y~x,se=T, level=0.95)+
        ggtitle(paste("r=",round(cor.test(orthVIP.new.LVD.adj, predVIP.new.LVD.adj, method = "pearson")$estimate, digit=2)
                      ,"p=",formatC(cor.test(orthVIP.new.LVD.adj, predVIP.new.LVD.adj, method = "pearson")$p.value,format="e",digits=2)))+
        theme_classic()+
        geom_hline(yintercept = 1.5,linetype="dashed", size=1, color="grey")+
        geom_vline(xintercept = 1.5, linetype="dashed", color="grey", size=1)

dev.off()


##add fold change to metabolite list
getfoldchange=function(meta, diet) {
        findCON=which(as.character(unlist(diet))=="CON")
        findTRT=which(as.character(unlist(diet))!="CON")
        n=ncol(meta)
        fold.change=rep(NA,n)
        names(fold.change)=colnames(meta)
        for (i in 1:n){
                con=sum(meta[findCON,i])
                trt=sum(meta[findTRT,i])
                if (con>trt) {fold.change[i]=-con/trt} 
                else {fold.change[i]=trt/con}
        }
        fold.change
}

foldchange.LVD.adj=getfoldchange(meta.new.LVD.adj, diet.LVD)
foldchange.LVD=getfoldchange(meta.new.LVD, diet.LVD)

sum(abs(foldchange.LVD.adj)>1.5, na.rm = TRUE) #52 metabolites with adj.foldchange>1.5
VIP_Fold_adj=list(names(which(VIP.LVD.adj>1.5)), names(which(abs(foldchange.LVD.adj)>1.5)))
names(VIP_Fold_adj)=c("VIP1.5","Fold1.5")
library(gplots)
pdf("fdrD_w_int_adjVIP&Fold_wo_int.pdf")
venn(VIP_Fold_adj)
sig.d.lvd<-temp.name[(Fdr.final.diet.LVD<0.05)&!is.na((Fdr.final.diet.LVD))]
VIP_Fold_sigD_adj=list(names(which(VIP.LVD.adj>1.5)), names(which(abs(foldchange.LVD.adj)>1.5)),sig.d.lvd)
names(VIP_Fold_sigD_adj)=c("VIP1.5","Fold1.5","sig.d")
venn(VIP_Fold_sigD_adj)
dev.off()

install.packages("devtools")
library(devtools)
install_github("ManonMartin/MBXUCL")
#Perform K-fold cross-validation and calculate cross-validation root mean square error (cvRMSE)
#https://rdrr.io/github/ManonMartin/MBXUCL/man/cvOPLSDA.html

#Compile regression results (interaction model) with VIP from strain-adjusted oplsda and strain-adjusted metabolite values
LVD.pthw.Adjvip.fold.adjMeta<-cbind(Meta.pathway, VIPvector, foldchange.LVD.adj, foldchange.LVD, t(Fdr.final.LVD),lm.transformation.LVD,lm.method.LVD,t(meta.LVD.adj))

#Compile regression results (interaction model) with VIP from strain-adjusted oplsda and strain-adjusted metabolite values
LVD.pthw.Adjvip.fold.rawMeta<-cbind(Meta.pathway, VIPvector, foldchange.LVD.adj, foldchange.LVD, t(Fdr.final.LVD),lm.transformation.LVD,lm.method.LVD,t(meta.LVD))

write.table(LVD.pthw.Adjvip.fold.adjMeta,"Output//LVD_pthw_lm_adjVIP_adjMet.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LVD.pthw.Adjvip.fold.rawMeta,"Output//LVD_pthw_lm_adjVIP_rawMet.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(data.LVD,"Output//LVD_masterdata_scaled&imputed.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)






