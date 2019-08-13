#Regression analysis for LVD vs. CON and for LP vs. CON for 8 strains (removing 6a_CC026, 6b_CC006)

rm(list = ls())
setwd("C:\\Users\\UNC NRI\\Desktop\\Jing NRI\\Metabolomics\\G0_081418_8strains3diets_JX")
masterdata<-read.csv("CC Metabolomics Raw Data File for Analysis_ver62317_ME.csv", header=F)

####Clean up column names#############################################
temp<-masterdata[,-(1:28)]# have a temperary subset of just the metabolites
temp.name<-as.character(unlist(temp[1,]))
head(temp.name)
len=length(temp.name)
#add<-rep("m", len)
#name.rev<-cbind(add,temp.name)
#add.m<-function(old.name){
#        paste(old.name[1],old.name[2],sep=".")
#}
#new.name<-apply(name.rev,1,add.m)#this is the 652 new names

data<-masterdata[-(1:1),]
colnames(data)=c(as.character(unlist(masterdata[1,1:28])),temp.name)#this is the cleaned data
dim(data) #160 680
summary(factor(as.character(data$Diet)))#four diet groups
#CON          Low Protein Low Vitamin D       ME 1%SS 
#38                38            42            42
summary(factor(as.character(data$RIX)))
#1a  1b  2a  2b  3a  3b  6a  6b  9a  9b 
#21  19  16  15  16  17  11   8  20  17

#######################################################################
#######(1)Now make a subset to remove MS diet, 6a&6b###################
#######################################################################
data.8s3d<-data[(data$Diet!="ME 1%SS")&(data$RIX!="6a")&(data$RIX!="6b"),]
dim(data.8s3d)
#[1] 105 680
summary(factor(as.character(data.8s3d$Diet)))
# CON   Low Protein Low Vitamin D 
# 35            33            37
summary(factor(as.character(data.8s3d$RIX)))
#1a 1b 2a 2b 3a 3b 9a 9b 
#15 15 11 11 12 14 15 12 
diet3<-factor(as.character(unlist(data.8s3d$Diet)))
strain8<-factor(as.character(unlist(data.8s3d$RIX)))
meta<-sapply(data.8s3d[,29:680],function(x) as.numeric(as.character(x)))
#[1] 105 652 "meta" is the 652 untransformed, but scaled&imputed, data for 3diets and 6strains.
n.met=ncol(meta)
bad.8s3d=apply(meta,2, function(x) (length(unique(x)) ==1 ))
sum(bad.8s3d)#4 metabolites with the same value, including NA, for all samples
which(bad.8s3d)
#m.4-vinylphenol sulfate              m.daidzein             m.genistein             m.glycitein 
#165                     268                     318                     339 

write.table(data.8s3d,"InputFile/imputedData_8strain3diet.txt", row.names = FALSE, quote=FALSE, sep = "\t")
#######################################################################
#######(2)Make 2 subsets, one for CON&LVD, the other for CON&LP########
#######################################################################
data.LVD<-data.8s3d[(data.8s3d$Diet!="Low Protein"),]
diet.LVD<-factor(as.character(unlist(data.LVD$Diet)))
strain.LVD<-factor(as.character(unlist(data.LVD$RIX)))
dim(data.LVD) #[1] 72 680
summary(data.LVD$Diet)
#CON          Diet   Low Protein Low Vitamin D       ME 1%SS 
#35             0             0            37             0 
summary(data.LVD$RIX)
#1a  1b  2a  2b  3a  3b  6a  6b  9a  9b RIX 
#10  10   8   7   8  10   0   0  10   9   0
meta.LVD<-sapply(data.LVD[,29:680],function(x) as.numeric(as.character(x)))
bad.LVD=apply(meta.LVD,2, function(x) (length(unique(x)) ==1 ))


data.LP<-data.8s3d[(data.8s3d$Diet!="Low Vitamin D"),]
diet.LP<-factor(as.character(unlist(data.LP$Diet)))
strain.LP<-factor(as.character(unlist(data.LP$RIX)))
dim(data.LP) #[1] 68 680
summary(data.LP$Diet)
#CON          Diet   Low Protein Low Vitamin D       ME 1%SS 
#35             0            33             0             0  
summary(data.LP$RIX)
#1a  1b  2a  2b  3a  3b  6a  6b  9a  9b RIX 
#10  10   7   8   7   8   0   0  10   8   0
meta.LP<-sapply(data.LP[,29:680],function(x) as.numeric(as.character(x)))
bad.LP=apply(meta.LP,2, function(x) (length(unique(x)) ==1 ))

write.table(data.LVD,"InputFile/imputedData_LVD.txt", row.names = FALSE, quote=FALSE, sep = "\t")
write.table(data.LP,"InputFile/imputedData_LP.txt", row.names = FALSE, quote=FALSE, sep = "\t")
#######################################################################
#######(3)CON&LVD: transformation to meet regression assumptions#######
#######################################################################
n.LVD=nrow(data.LVD) #number of samples in CON&LVD
library(car)
options(contrasts = c("contr.sum", "contr.poly")) 
library(lmtest)# bptest for testing heteroscedasticity


######(3-a) no transformation########
res.meta.LVD<-matrix(data=NA,nrow=n.LVD, ncol=n.met)
colnames(res.meta.LVD)=temp.name
bp.meta.LVD<-rep(NA,length=n.met)
sw.meta.LVD<-rep(NA,length=n.met)
pb <- winProgressBar(title="R Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i]) next
        res.meta.LVD[,i]=resid(lm(meta.LVD[,i]~diet.LVD*strain.LVD))
        bp.meta.LVD[i]=bptest(lm(meta.LVD[,i]~diet.LVD*strain.LVD))$p.value
        sw.meta.LVD[i]=shapiro.test(res.meta.LVD[,i])$p.value
        
        info<-sprintf("%d%% done", round(i/n.met*100))
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum(bp.meta.LVD>=0.05,na.rm=TRUE) #391 homostasciticity 
sum(bp.meta.LVD<0.05,na.rm=TRUE) #257 heterostasciticity 
sum(sw.meta.LVD>=0.05,na.rm=TRUE) #174 normal 
sum(sw.meta.LVD<0.05,na.rm=TRUE) #474 not normally distributed 

sum((bp.meta.LVD>=0.05)&(sw.meta.LVD>=0.05),na.rm=TRUE) #107 meet both requirements without transformation

######(3-b) BoxCox transformation########
library(MASS)
meta.box.LVD=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.box.LVD)=temp.name
#http://rcompanion.org/handbook/I_12.html# 
corticosterone.LVD.raw=lm(meta.LVD[,which(temp.name=="corticosterone")]~diet.LVD*strain.LVD)
Box=boxcox(corticosterone.LVD.raw,lambda = seq(-6,6,0.1), plotit = TRUE)
Cox=data.frame(Box$x,Box$y) 
Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
lambda= Cox2[1, "Box.x"]# Extract the lambda with greatest log-likelihood  

pb <- winProgressBar(title="Transformation Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i])next
        model=lm(meta.LVD[,i]~diet.LVD*strain.LVD)
        Box=boxcox(model,lambda = seq(-6,6,0.1))
        Cox=data.frame(Box$x,Box$y) 
        Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
        lambda= Cox2[1, "Box.x"]# Extract that lambda
        if(lambda==0){meta.box.LVD[,i]=log(meta.LVD[,i])}else{
                meta.box.LVD[,i]=(meta.LVD[,i] ^ lambda - 1)/lambda
        }
        
        rm(Box,Cox,Cox2,lambda)# or rm(list=c('Box','Cox'))
        
        info<-sprintf("%d th metabolite transformed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

res.box.LVD<-matrix(data=NA,nrow= n.LVD, ncol=n.met)
colnames(res.box.LVD)=temp.name
bp.box.LVD<-rep(NA,length=n.met)
sw.box.LVD<-rep(NA,length=n.met)

pb <- winProgressBar(title="Regression Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i])next
        res.box.LVD[,i]=resid(lm(meta.box.LVD[,i]~diet.LVD*strain.LVD))
        bp.box.LVD[i]=bptest(lm(meta.box.LVD[,i]~diet.LVD*strain.LVD))$p.value
        sw.box.LVD[i]=shapiro.test(res.box.LVD[,i])$p.value

        info<-sprintf("%d th metabolite done", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum((bp.box.LVD>=0.05)&(sw.box.LVD>=0.05),na.rm=TRUE) #412/652 meet both requirements after BoxCox
sum(((bp.box.LVD>=0.05)&(sw.box.LVD>=0.05))|((bp.meta.LVD>=0.05)&(sw.meta.LVD>=0.05))
    ,na.rm=TRUE) #425/652 meet both requirements before or after BoxCox

######(3-c) Tukey's Ladder of Powers transformation########
meta.tuk.LVD=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.tuk.LVD)=temp.name
res.tuk.LVD<-matrix(data=NA,nrow= n.LVD, ncol=n.met)
colnames(res.tuk.LVD)=temp.name
bp.tuk.LVD<-rep(NA,length=n.met)
sw.tuk.LVD<-rep(NA,length=n.met)

library(rcompanion)#http://rcompanion.org/handbook/I_12.html 
pb <- winProgressBar(title="Transformation Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LVD[i])next
        meta.tuk.LVD[,i]=transformTukey(meta.LVD[,i],plotit=FALSE,verbose = FALSE)
        
        info<-sprintf("%d th metabolite transformed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

for(i in 1:n.met){
        if(bad.LVD[i]) next
        res.tuk.LVD[,i]=resid(lm(meta.tuk.LVD[,i]~diet.LVD*strain.LVD))
        bp.tuk.LVD[i]=bptest(lm(meta.tuk.LVD[,i]~diet.LVD*strain.LVD))$p.value
        sw.tuk.LVD[i]=shapiro.test(res.tuk.LVD[,i])$p.value
}
sum((bp.tuk.LVD>=0.05)&(sw.tuk.LVD>=0.05),na.rm=TRUE) #401/652 meet both requirements after Tukey transformation
sum(((bp.box.LVD>=0.05)&(sw.box.LVD>=0.05))|((bp.meta.LVD>=0.05)&(sw.meta.LVD>=0.05))
    |((bp.tuk.LVD>=0.05)&(sw.tuk.LVD>=0.05))
    ,na.rm=TRUE) #445/652 meet both requirements before or after BoxCox or Tukey transformation

######(3-d) Log transformation########
meta.log.LVD=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.log.LVD)=temp.name
res.log.LVD<-matrix(data=NA,nrow= n.LVD, ncol=n.met)
colnames(res.log.LVD)=temp.name
bp.log.LVD<-rep(NA,length=n.met)
sw.log.LVD<-rep(NA,length=n.met)

for (i in 1:n.met){
        if(bad.LVD[i])next
        meta.log.LVD[,i]=log(meta.LVD[,i])
        res.log.LVD[,i]=resid(lm(meta.log.LVD[,i]~diet.LVD*strain.LVD))
        bp.log.LVD[i]=bptest(lm(meta.log.LVD[,i]~diet.LVD*strain.LVD))$p.value
        sw.log.LVD[i]=shapiro.test(res.log.LVD[,i])$p.value
}
sum((bp.log.LVD>=0.05)&(sw.log.LVD>=0.05),na.rm=TRUE) #307/652 meet both requirements after Tukey transformation
sum(((bp.box.LVD>=0.05)&(sw.box.LVD>=0.05))|((bp.meta.LVD>=0.05)&(sw.meta.LVD>=0.05))
    |((bp.tuk.LVD>=0.05)&(sw.tuk.LVD>=0.05))|((bp.log.LVD>=0.05)&(sw.log.LVD>=0.05))
    ,na.rm=TRUE) #454/652 meet both requirements before or after BoxCox or Tukey or log transformation

######(3-e) combine transformation files########
Pval_assump_LVD<-rbind(bp.meta.LVD,sw.meta.LVD,bp.box.LVD,sw.box.LVD,bp.tuk.LVD,sw.tuk.LVD,bp.log.LVD,sw.log.LVD)
colnames(Pval_assump_LVD)=temp.name
rawGood.LVD<-apply(Pval_assump_LVD[1:2,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
boxGood.LVD<-apply(Pval_assump_LVD[3:4,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
tukGood.LVD<-apply(Pval_assump_LVD[5:6,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
logGood.LVD<-apply(Pval_assump_LVD[7:8,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))

apply(Pval_assump_LVD,1,function(x) sum(is.na(x)))
#bp.meta.LVD sw.meta.LVD  bp.box.LVD  sw.box.LVD  bp.tuk.LVD  sw.tuk.LVD  bp.log.LVD  sw.log.LVD 
#4           4           4           4           4           4           4           4 

Good.LVD=rawGood.LVD|logGood.LVD|boxGood.LVD|tukGood.LVD
sum(Good.LVD,na.rm = TRUE) #454
sum(!Good.LVD,na.rm = TRUE) #194
sum(is.na(Good.LVD)) #4
Compliance.LVD<-rbind(Good.LVD,rawGood.LVD, boxGood.LVD, tukGood.LVD,logGood.LVD)
colnames(Compliance.LVD)=temp.name

#priority high-low: raw, boxcox, tukey, log
data.choice.LVD=vector(length=n.met)
for (i in 1:n.met){
        if(bad.LVD[i]){data.choice.LVD[i]="noChoice"}
        else if(rawGood.LVD[i]==TRUE) {data.choice.LVD[i]="raw"}
        else if(boxGood.LVD[i]==TRUE) {data.choice.LVD[i]="boxcox"}
        else if(tukGood.LVD[i]==TRUE) {data.choice.LVD[i]="tuk"}
        else if(logGood.LVD[i]==TRUE) {data.choice.LVD[i]="log"}
        else {data.choice.LVD[i]="noChoice"}
}
table(data.choice.LVD)
#boxcox      log noChoice      raw      tuk 
#318        9      198      107       20 
bp=barplot(table(data.choice.LVD),main="transformation for OLS-linear regression", 
           xlab="type of transformation", ylab="number of metabolites")
text(bp,15,table(data.choice.LVD), cex=1, pos=3)

bh.LVD<-apply(Pval_assump_LVD[c(1,3,5,7),],2, function(x) (x[1]>=0.05)|(x[2]>=0.05)|(x[3]>=0.05)|(x[4]>=0.05))
sw.LVD<-apply(Pval_assump_LVD[c(2,4,6,8),],2, function(x) (x[1]>=0.05)|(x[2]>=0.05)|(x[3]>=0.05)|(x[4]>=0.05))

sum((data.choice.LVD=="noChoice")&!sw.LVD&!bh.LVD,na.rm=TRUE)#17 meet neither, plus 4 bad ==>21
sum((data.choice.LVD=="noChoice")&sw.LVD&bh.LVD,na.rm=TRUE)#22 meet both, but under different transformation

sum((data.choice.LVD=="noChoice")&sw.LVD,na.rm=TRUE)#107/198 normal but heteroscedastic, use robust errors (HC errors)
sum((data.choice.LVD=="noChoice")&bh.LVD,na.rm=TRUE)#92/198 homoscedastic but non-normal
sum((data.choice.LVD=="noChoice")&bh.LVD&!sw.LVD,na.rm=TRUE)#70/198 homoscedastic but non-normal, use robust regression
#In conclusion, among the 198 no-choice, 107 with OLS-robust-error, the rest 70 use robust regression, the rest 21 toss.
#OLS-regression for all 454 metabolites that meet both assumptions w/ or w/o transformation.


#Now specify data transformation to use for the 198.The 4 NA metabolites are included but will be properly excluded 
noChoice.HC.LVD=(data.choice.LVD=="noChoice"&sw.LVD) #specify positions for 107 metabolites for HC/robust errors
noChoice.rr.LVD=(data.choice.LVD=="noChoice"&bh.LVD&!sw.LVD) #specify positions for the 70 metabolites for robust regression
sum(noChoice.rr.LVD&noChoice.HC.LVD,na.rm=TRUE) #0, confirm the two don't overlap
#determine which transformation to do, priority high-low: raw, boxcox, tukey, log
#Pval_assump_LVD<-rbind(bp.meta.LVD,sw.meta.LVD,bp.box.LVD,sw.box.LVD,bp.tuk.LVD,sw.tuk.LVD,bp.log.LVD,sw.log.LVD)
data.choice.HC.LVD=rep(NA,length=n.met)
for (i in 1:n.met){
        if (!noChoice.HC.LVD[i]|bad.LVD[i]) next
        if(Pval_assump_LVD[2,i]>=0.05) {data.choice.HC.LVD[i]="raw"}
        else if(Pval_assump_LVD[4,i]>=0.05) {data.choice.HC.LVD[i]="boxcox"}
        else if(Pval_assump_LVD[6,i]>=0.05) {data.choice.HC.LVD[i]="tuk"}
        else if(Pval_assump_LVD[8,i]>=0.05) {data.choice.HC.LVD[i]="log"}
}
table(data.choice.HC.LVD)
#boxcox    log    raw    tuk 
# 67      1     36      3 


data.choice.rr.LVD=rep(NA,length=n.met)
for (i in 1:n.met){
        if ((!noChoice.rr.LVD[i])|bad.LVD[i]) next
        if(Pval_assump_LVD[1,i]>=0.05) {data.choice.rr.LVD[i]="raw"}
        else if(Pval_assump_LVD[3,i]>=0.05) {data.choice.rr.LVD[i]="boxcox"}
        else if(Pval_assump_LVD[5,i]>=0.05) {data.choice.rr.LVD[i]="tuk"}
        else if(Pval_assump_LVD[7,i]>=0.05) {data.choice.rr.LVD[i]="log"}
}
table(data.choice.rr.LVD)
#boxcox    log    raw 
#16      5     49 
sum(!is.na(data.choice.rr.LVD))#70

#######################################################################
#######(4)CON&LVD: regression analysis#################################
#######################################################################
meta.mix.LVD=matrix(nrow=n.LVD,ncol=n.met)
colnames(meta.mix.LVD)=temp.name
dim(meta.mix.LVD) #[1]  72 652
meta.mix.LVD[,data.choice.LVD=="raw"]=meta.LVD[,data.choice.LVD=="raw"]#cannot use this when vector has na values
meta.mix.LVD[,data.choice.LVD=="boxcox"]=meta.box.LVD[,data.choice.LVD=="boxcox"]
meta.mix.LVD[,data.choice.LVD=="tuk"]=meta.tuk.LVD[,data.choice.LVD=="tuk"]
meta.mix.LVD[,data.choice.LVD=="log"]=meta.log.LVD[,data.choice.LVD=="log"]

meta.mix.LVD[,which(data.choice.HC.LVD=="raw")]=meta.LVD[,which(data.choice.HC.LVD=="raw")]#because which() can filter na
meta.mix.LVD[,which(data.choice.HC.LVD=="boxcox")]=meta.box.LVD[,which(data.choice.HC.LVD=="boxcox")]
meta.mix.LVD[,which(data.choice.HC.LVD=="tuk")]=meta.tuk.LVD[,which(data.choice.HC.LVD=="tuk")]
meta.mix.LVD[,which(data.choice.HC.LVD=="log")]=meta.log.LVD[,which(data.choice.HC.LVD=="log")]

meta.mix.LVD[,which(data.choice.rr.LVD=="raw")]=meta.LVD[,which(data.choice.rr.LVD=="raw")]
meta.mix.LVD[,which(data.choice.rr.LVD=="boxcox")]=meta.box.LVD[,which(data.choice.rr.LVD=="boxcox")]
meta.mix.LVD[,which(data.choice.rr.LVD=="tuk")]=meta.tuk.LVD[,which(data.choice.rr.LVD=="tuk")]
meta.mix.LVD[,which(data.choice.rr.LVD=="log")]=meta.log.LVD[,which(data.choice.rr.LVD=="log")]

sum(apply(meta.mix.LVD,2,function(x) sum(is.na(x))))/n.LVD #Confirmed 21 metabolites not added to data for further analysis

#######################################################################
#######(4-a)CON&LVD: OLS linear regression for 454 metabolites#########
#######################################################################
library(car)
options(contrasts = c("contr.sum", "contr.poly")) 
test.model=lm(meta.mix.LVD[,1]~diet.LVD*strain.LVD)
summary(test.model)$adj.r.squared #[1] 0.6167251
#http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
#http://blog.minitab.com/blog/adventures-in-statistics-2/what-is-the-f-test-of-overall-significance-in-regression-analysis 
lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
}
lmp(test.model)#[1] 9.205171e-10, this is the pval on the F-test of lm 

test=Anova(test.model, type="III")
test.p<-test[,4]
names(test.p)=row.names(test)# extract pval from ANOVA with names.

OLS.p.LVD<-matrix(nrow=length(test.p),ncol=n.met)
colnames(OLS.p.LVD)=temp.name
row.names(OLS.p.LVD)=row.names(test)

pb <- winProgressBar(title="Anova progress", min=0, max=100, initial=0,label = paste("0% done"))
for(i in 1:n.met){
        if (data.choice.LVD[i]=="noChoice") next
        model=lm(meta.mix.LVD[,i]~diet.LVD*strain.LVD)
        test=Anova(model, type="III")
        OLS.p.LVD[,i]=test[,4]
        rm(model, test)
        
        info<-sprintf("%d th metabolite analyzed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}
#Last 4 metabolites are noChoice, so don't freak out if run stops at 648th metabolite.

#######################################################################
#######(4-b)CON&LVD: heteroskedasticity-consistent (HC) errors for 107 metabolites
#######################################################################
#Robust standard errors following lm() can deal with the 107 metabolites in meta.mix.LVD[,noChoice.HC.LVD] 
#Calculate ANOVA significance and add to the anova.p dataframe.

library(car)
options(contrasts = c("contr.sum", "contr.poly")) 
#http://thestatsgeek.com/2014/02/14/the-robust-sandwich-variance-estimator-for-linear-regression-using-r/
#https://stackoverflow.com/questions/4385436/regression-with-heteroskedasticity-corrected-standard-errors 
#https://stats.stackexchange.com/questions/131401/how-to-get-anova-table-with-robust-standard-errors/132521 ##Wald for robust errors
#http://data.princeton.edu/wws509/notes/c2s3.html ##Basis of Wald test being identical to ANOVA for linear model
#http://data.princeton.edu/wws509/r/robust.html ##Robust standard errors in R
library(sandwich)
library(lmtest)

corticosterone.LVD=lm(meta.mix.LVD[,which(temp.name=="corticosterone")]~diet.LVD*strain.LVD)
Anova(corticosterone.LVD,type=3)# OLS without HC standard error
data.choice.HC.LVD[which(temp.name=="corticosterone")]#BoxCox

sandwich_se.LVD <- diag(vcovHC(corticosterone.LVD, type = "HC1"))^0.5
sandwich_se.LVD  ##this is the robust standard errors, same with those in coeftest(, vcov.=vcovHC) ##
coeftest(corticosterone.LVD, vcov. = vcovHC(corticosterone.LVD,type = "HC1"))
Anova(corticosterone.LVD,type=3,vcov. = vcovHC(corticosterone.LVD,type = "HC1"))$`Pr(>F)`
#Coefficient covariances computed by vcovHC(corticosterone.LVD, type = "HC1")
#[1] 1.048401e-05 2.909330e-01 3.373502e-06 7.409205e-01           NA
#[1] "(Intercept)"         "diet.LVD"            "strain.LVD"          "diet.LVD:strain.LVD" "Residuals"

row.names(OLS.p.LVD) #[1] "(Intercept)"         "diet.LVD"            "strain.LVD"          "diet.LVD:strain.LVD" "Residuals"
pb <- winProgressBar(title="Wald with robust error progress", min=0, max=100, initial=0,label = paste("0% done"))
for(i in 1:n.met){
        if ((!noChoice.HC.LVD[i])|is.na(noChoice.HC.LVD[i])) next
        model=lm(meta.mix.LVD[,i]~diet.LVD*strain.LVD)
        OLS.p.LVD[,i]=Anova(model,type=3,vcov. = vcovHC(model,type = "HC1"))$`Pr(>F)`
        rm(model)
        
        info<-sprintf("%d th metabolite analyzed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}
OLS.p.LVD[,which(temp.name=="corticosterone")] ##to confirm loop runned correctly
# (Intercept)            diet.LVD          strain.LVD diet.LVD:strain.LVD           Residuals 
#1.048401e-05        2.909330e-01        3.373502e-06        7.409205e-01                  NA 

#######################################################################
#######(4-c)CON&LVD: robust regression for 70 metabolites
#######################################################################
library(MASS)
#Mass:rml has issue with list that has unused levels. Convert list to character.
library(car)
options(contrasts = c("contr.sum", "contr.poly")) 

rr.Smethylmethionine <- rlm(meta.mix.LVD[,which(temp.name=="S-methylmethionine")]~diet.LVD*strain.LVD)
summary(rr.Smethylmethionine)
x=linearHypothesis(rr.Smethylmethionine, c("diet.LVD1")) 
y=linearHypothesis(rr.Smethylmethionine, c("strain.LVD1","strain.LVD2","strain.LVD3","strain.LVD4","strain.LVD5","strain.LVD6"
                                        ,"strain.LVD7")) 
z=linearHypothesis(rr.Smethylmethionine, c("diet.LVD1:strain.LVD1","diet.LVD1:strain.LVD2"
                                        ,"diet.LVD1:strain.LVD3","diet.LVD1:strain.LVD4"
                                        ,"diet.LVD1:strain.LVD5","diet.LVD1:strain.LVD6"
                                        ,"diet.LVD1:strain.LVD7"))
rr.Smethylmethionine.p=c(x$`Pr(>F)`[2],y$`Pr(>F)`[2],z$`Pr(>F)`[2])
#[1] 8.745488e-163 1.515547e-198 1.191488e-170

test1=lm(meta.mix.LVD[,which(temp.name=="S-methylmethionine")]~diet.LVD*strain.LVD)
Anova(test1,type=3)$`Pr(>F)`        
#[1] 4.188756e-25 5.160799e-02 4.420928e-09 6.625636e-01           NA

#Start a different pvalue table for comparing the two versions for the 70.
sensit.ols.LVD=OLS.p.LVD
sensit.rr.LVD=OLS.p.LVD
row.names(sensit.rr.LVD)
#[1] "(Intercept)"         "diet.LVD"            "strain.LVD"          "diet.LVD:strain.LVD" "Residuals"    

pb <- winProgressBar(title="Sensitivity of the 70 metabolites", min=0, max=100, initial=0,label = paste("0% done"))
for(i in 1:n.met){
        if ((!noChoice.rr.LVD[i])|is.na(noChoice.rr.LVD[i])) next
        rr.model=rlm(meta.mix.LVD[,i]~diet.LVD*strain.LVD,maxit = 200)
        sensit.rr.LVD[2,i]=linearHypothesis(rr.model, c("diet.LVD1"))$`Pr(>F)`[2] 
        sensit.rr.LVD[3,i]=linearHypothesis(rr.model, c("strain.LVD1","strain.LVD2","strain.LVD3"
                                                       ,"strain.LVD4","strain.LVD5","strain.LVD6"
                                                       ,"strain.LVD7"))$`Pr(>F)`[2] 
        sensit.rr.LVD[4,i]=linearHypothesis(rr.model, c("diet.LVD1:strain.LVD1","diet.LVD1:strain.LVD2"
                                                       ,"diet.LVD1:strain.LVD3","diet.LVD1:strain.LVD4"
                                                       ,"diet.LVD1:strain.LVD5","diet.LVD1:strain.LVD6"
                                                       ,"diet.LVD1:strain.LVD7"))$`Pr(>F)`[2]
        ols.model=lm(meta.mix.LVD[,i]~diet.LVD*strain.LVD)
        test=Anova(ols.model, type="III")
        sensit.ols.LVD[,i]=test[,4]
        
        rm(rr.model,ols.model, test)
        
        info<-sprintf("%d th metabolite analyzed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum(is.na(sensit.ols.LVD[2,which(noChoice.rr.LVD)])) #0
sum(is.na(sensit.rr.LVD[2,which(noChoice.rr.LVD)])) #0
sensit.rr.LVD[,which(temp.name=="S-methylmethionine")]
#(Intercept)            diet.LVD          strain.LVD diet.LVD:strain.LVD           Residuals 
#NA       8.745488e-163       1.515547e-198       1.191488e-170                  NA 

sensi.diet.LVD<-rep(NA,n.met)
sensi.strain.LVD<-rep(NA,n.met)
sensi.int.LVD<-rep(NA,n.met)
sensi.LVD=rbind(sensi.diet.LVD,sensi.strain.LVD,sensi.strain.LVD)
row.names(sensi.LVD)=c("diet","strain","diet:strain")
colnames(sensi.LVD)=temp.name

##Filter robust regression results so that if significant only present in rr, not in OLS, the result is removed.
for(i in 1:n.met){
        if ((!noChoice.rr.LVD[i])|is.na(noChoice.rr.LVD[i])) next
        
        if((sensit.ols.LVD[2,i]>0.05)&(sensit.rr.LVD[2,i]<0.05)){sensi.LVD[1,i]="Exclude"} 
        else{sensi.LVD[1,i]="Pass"}
        if((sensit.ols.LVD[3,i]>0.05)&(sensit.rr.LVD[3,i]<0.05)){sensi.LVD[2,i]="Exclude"} 
        else{sensi.LVD[2,i]="Pass"}
        if((sensit.ols.LVD[4,i]>0.05)&(sensit.rr.LVD[4,i]<0.05)){sensi.LVD[3,i]="Exclude"} 
        else{sensi.LVD[3,i]="Pass"}
}

sens.3effects.LVD<-apply(sensi.LVD,2,function(x) x[1]=="Pass"&x[2]=="Pass"&x[3]=="Pass")
sum(sens.3effects.LVD,na.rm=TRUE)
#[1] 31 metabolites passed sensitivity test and the robust regression results for those were kept

#These are metabolites passing sensitivity test with significant diet effect from Robust Regression
pass.sens.diet.rr.LVD<-temp.name[which((sens.3effects.LVD==TRUE)&sensit.rr.LVD[2,]<0.05)]
#[1] "1-oleoyl-2-linoleoyl-GPG (18:1/18:2)*" "1-oleoyl-GPG (18:1)*"                  "cytidine" 

##Merge p-values from Robust Regression that passes sensitivity test to the master p-value table
Pval.final.LVD=OLS.p.LVD[2:4,]
row.names(Pval.final.LVD)=c("diet.pval","strain.pval", "interaction.pval")
for (i in 1:n.met) {
        if (sens.3effects.LVD[i]&!is.na(sens.3effects.LVD[i])){
                Pval.final.LVD[,i]=sensit.rr.LVD[2:4,i]
        }
}

#######################################################################
#######(4-D)CON&LVD: Compile regression results########################
#######################################################################
Pval.final.LVD[,which(sens.3effects.LVD==FALSE)]##Make sure those failed sensitivity all had NA Pval
sum(is.na(Pval.final.LVD[1,]))#[1] 60 not analyzed, including 21 not meeting assumptions and 39 failed sensitivity test

Fdr.final.diet.LVD=p.adjust(Pval.final.LVD[1,],method="fdr", n=sum(!is.na(Pval.final.LVD[1,])))#FDR with 592 metabolites
Fdr.final.strain.LVD=p.adjust(Pval.final.LVD[2,],method="fdr", n=sum(!is.na(Pval.final.LVD[2,])))
Fdr.final.int.LVD=p.adjust(Pval.final.LVD[3,],method="fdr", n=sum(!is.na(Pval.final.LVD[3,])))
Fdr.final.LVD=rbind(Fdr.final.diet.LVD,Fdr.final.strain.LVD,Fdr.final.int.LVD)
row.names(Fdr.final.LVD)=c("diet.fdr","strain.fdr", "interaction.fdr")

colnames(Fdr.final.LVD)=temp.name
lm.method.LVD=vector(mode="character", length=n.met)
lm.method.LVD[which(data.choice.LVD!="noChoice")]="OLS-assumptionsMet"
lm.method.LVD[which(noChoice.HC.LVD)]="HC-error"
lm.method.LVD[which(sens.3effects.LVD==TRUE)]="RobustRegressionPassingSensitivity"
lm.method.LVD[which(sens.3effects.LVD==FALSE)]="FailingSensitivity"

table(lm.method.LVD)
#                        FailingSensitivity                           HC-error 
#21                                 39                                107 
#OLS-assumptionsMet RobustRegressionPassingSensitivity 
#454                                 31
op <- par(mar=c(15,4,4,2))
methodlist=table(lm.method.LVD)[c(4,3,5,2,1)]
barname=names(methodlist)
barname[which(barname=="")]="NotAnalyzed"
bp=barplot(methodlist,main="Regression methods", ylim=c(0,max(methodlist)+100),
           ylab="number of metabolites", names.arg = barname, horiz = F, las=2)
text(bp,methodlist,methodlist, cex=1, pos=3)

lm.transformation.LVD=vector(mode="character", length=n.met)
lm.transformation.LVD[which(data.choice.LVD!="noChoice")]=data.choice.LVD[which(data.choice.LVD!="noChoice")]
table(lm.transformation.LVD)
#     boxcox    log    raw    tuk 
#198    318      9    107     20
lm.transformation.LVD[!is.na(data.choice.rr.LVD)]=data.choice.rr.LVD[!is.na(data.choice.rr.LVD)]
table(lm.transformation.LVD)
#       boxcox    log    raw    tuk 
#128    334     14    156     20
lm.transformation.LVD[!is.na(data.choice.HC.LVD)]=data.choice.HC.LVD[!is.na(data.choice.HC.LVD)]
table(lm.transformation.LVD)
#       boxcox    log    raw    tuk 
#21    401     15    192     23

op <- par(mar=c(15,4,4,2))
transformation=table(lm.transformation.LVD)[c(4,2,5,3,1)]
barname=names(transformation)
barname[which(barname=="")]="NotAnalyzed"
bp=barplot(transformation,main="Transformation methods", ylim=c(0,max(transformation)+100),
           ylab="number of metabolites", names.arg = barname, horiz = F, las=2)
text(bp,transformation,transformation, cex=1, pos=3)
rm(op)

##The compiled results table has 652 columns with each col corresponding to a metabolite##
##Row1-3 have p-values for diet, strain and diet*strain interaction##
##Row4-6 have fdr for diet, strain and diet*strain interaction##
##Row7 tells if data was transformed (raw for no transformation) and which transformation is used (box-cox, tukey or log)##
##Row8 tells which regression method was used:
#### 454 metabolites met both residual normality and residual homoscedasticity, therefore ordinary least squares (OLS) regression and anova to test effects.
#### 107 metabolites only met residual normality, therefore linear regression lwith heteroskedasticity-consistent standard errors and Wald test to test effects.
#### 70 metabolites only met residual homoscedasticity, therefore robust regression and exclusion of 39 metabolites that failed sensitivity test (significance in robust regression but not in OLS)

CompiledResults_LVD<-rbind(Pval.final.LVD,Fdr.final.LVD,lm.transformation.LVD,lm.method.LVD)
colnames(CompiledResults_LVD)=temp.name
#This is the metabolite data table with proper transformation for the regression analysis
Data.wTransformation.LVD<-rbind(meta.mix.LVD, CompiledResults_LVD)
Data.woTransformation.LVD<-rbind(meta.LVD, CompiledResults_LVD)
#Animal sample information is "CC Metabolomics Raw Data File for Analysis_ver62317_ME".

##None of the output files have col names. Use col names of "CC Metabolomics Raw Data File for Analysis_ver62317_ME". 
write.table(Data.wTransformation.LVD,"LVD_outputFile/LVD_regressionResult_dataTransformed.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(Data.woTransformation.LVD,"LVD_outputFile/LVD_regressionResult_dataNotTransformed.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)








#######################################################################
#######(5)CON&LP: transformation to meet regression assumptions#######
#######################################################################
n.LP=nrow(data.LP) #number of samples in CON&LP
library(car)
options(contrasts = c("contr.sum", "contr.poly")) 
library(lmtest)# bptest for testing heteroscedasticity


######(5-a) no transformation########
res.meta.LP<-matrix(data=NA,nrow=n.LP, ncol=n.met)
colnames(res.meta.LP)=temp.name
bp.meta.LP<-rep(NA,length=n.met)
sw.meta.LP<-rep(NA,length=n.met)
pb <- winProgressBar(title="R Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LP[i]) next
        res.meta.LP[,i]=resid(lm(meta.LP[,i]~diet.LP*strain.LP))
        bp.meta.LP[i]=bptest(lm(meta.LP[,i]~diet.LP*strain.LP))$p.value
        sw.meta.LP[i]=shapiro.test(res.meta.LP[,i])$p.value
        
        info<-sprintf("%d%% done", round(i/n.met*100))
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum(bp.meta.LP>=0.05,na.rm=TRUE) #308 homostasciticity 
sum(bp.meta.LP<0.05,na.rm=TRUE) #340 heterostasciticity 
sum(sw.meta.LP>=0.05,na.rm=TRUE) #183 normal 
sum(sw.meta.LP<0.05,na.rm=TRUE) #465 not normally distributed 

sum((bp.meta.LP>=0.05)&(sw.meta.LP>=0.05),na.rm=TRUE) #100 meet both requirements without transformation

######(5-b) BoxCox transformation########
library(MASS)
meta.box.LP=matrix(nrow=n.LP,ncol=n.met)
colnames(meta.box.LP)=temp.name
pb <- winProgressBar(title="Transformation Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LP[i])next
        model=lm(meta.LP[,i]~diet.LP*strain.LP)
        Box=boxcox(model,lambda = seq(-6,6,0.1))
        Cox=data.frame(Box$x,Box$y) 
        Cox2 = Cox[with(Cox, order(-Cox$Box.y)),]
        lambda= Cox2[1, "Box.x"]# Extract that lambda
        if(lambda==0){meta.box.LP[,i]=log(meta.LP[,i])}else{
                meta.box.LP[,i]=(meta.LP[,i] ^ lambda - 1)/lambda
        }
        
        rm(Box,Cox,Cox2,lambda)# or rm(list=c('Box','Cox'))
        
        info<-sprintf("%d th metabolite transformed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

res.box.LP<-matrix(data=NA,nrow= n.LP, ncol=n.met)
colnames(res.box.LP)=temp.name
bp.box.LP<-rep(NA,length=n.met)
sw.box.LP<-rep(NA,length=n.met)

pb <- winProgressBar(title="Regression Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LP[i])next
        res.box.LP[,i]=resid(lm(meta.box.LP[,i]~diet.LP*strain.LP))
        bp.box.LP[i]=bptest(lm(meta.box.LP[,i]~diet.LP*strain.LP))$p.value
        sw.box.LP[i]=shapiro.test(res.box.LP[,i])$p.value
        
        info<-sprintf("%d th metabolite done", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum((bp.box.LP>=0.05)&(sw.box.LP>=0.05),na.rm=TRUE) #391/652 meet both requirements after BoxCox
sum(((bp.box.LP>=0.05)&(sw.box.LP>=0.05))|((bp.meta.LP>=0.05)&(sw.meta.LP>=0.05))
    ,na.rm=TRUE) #398/652 meet both requirements before or after BoxCox

######(5-c) Tukey's Ladder of Powers transformation########
meta.tuk.LP=matrix(nrow=n.LP,ncol=n.met)
colnames(meta.tuk.LP)=temp.name
res.tuk.LP<-matrix(data=NA,nrow= n.LP, ncol=n.met)
colnames(res.tuk.LP)=temp.name
bp.tuk.LP<-rep(NA,length=n.met)
sw.tuk.LP<-rep(NA,length=n.met)

library(rcompanion)#http://rcompanion.org/handbook/I_12.html 
pb <- winProgressBar(title="Transformation Progress", min=0, max=100, initial=0,label = paste("0% done"))
for (i in 1:n.met){
        if(bad.LP[i])next
        meta.tuk.LP[,i]=transformTukey(meta.LP[,i],plotit=FALSE,verbose = FALSE)
        
        info<-sprintf("%d th metabolite transformed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

for(i in 1:n.met){
        if(bad.LP[i]) next
        res.tuk.LP[,i]=resid(lm(meta.tuk.LP[,i]~diet.LP*strain.LP))
        bp.tuk.LP[i]=bptest(lm(meta.tuk.LP[,i]~diet.LP*strain.LP))$p.value
        sw.tuk.LP[i]=shapiro.test(res.tuk.LP[,i])$p.value
}
sum((bp.tuk.LP>=0.05)&(sw.tuk.LP>=0.05),na.rm=TRUE) #369/652 meet both requirements after Tukey transformation
sum(((bp.box.LP>=0.05)&(sw.box.LP>=0.05))|((bp.meta.LP>=0.05)&(sw.meta.LP>=0.05))
    |((bp.tuk.LP>=0.05)&(sw.tuk.LP>=0.05))
    ,na.rm=TRUE) #412/652 meet both requirements before or after BoxCox or Tukey transformation

######(5-d) Log transformation########
meta.log.LP=matrix(nrow=n.LP,ncol=n.met)
colnames(meta.log.LP)=temp.name
res.log.LP<-matrix(data=NA,nrow= n.LP, ncol=n.met)
colnames(res.log.LP)=temp.name
bp.log.LP<-rep(NA,length=n.met)
sw.log.LP<-rep(NA,length=n.met)

for (i in 1:n.met){
        if(bad.LP[i])next
        meta.log.LP[,i]=log(meta.LP[,i])
        res.log.LP[,i]=resid(lm(meta.log.LP[,i]~diet.LP*strain.LP))
        bp.log.LP[i]=bptest(lm(meta.log.LP[,i]~diet.LP*strain.LP))$p.value
        sw.log.LP[i]=shapiro.test(res.log.LP[,i])$p.value
}
sum((bp.log.LP>=0.05)&(sw.log.LP>=0.05),na.rm=TRUE) #282/652 meet both requirements after Tukey transformation
sum(((bp.box.LP>=0.05)&(sw.box.LP>=0.05))|((bp.meta.LP>=0.05)&(sw.meta.LP>=0.05))
    |((bp.tuk.LP>=0.05)&(sw.tuk.LP>=0.05))|((bp.log.LP>=0.05)&(sw.log.LP>=0.05))
    ,na.rm=TRUE) #420/652 meet both requirements before or after BoxCox or Tukey or log transformation

######(5-e) combine transformation files########
Pval_assump_LP<-rbind(bp.meta.LP,sw.meta.LP,bp.box.LP,sw.box.LP,bp.tuk.LP,sw.tuk.LP,bp.log.LP,sw.log.LP)
colnames(Pval_assump_LP)=temp.name
rawGood.LP<-apply(Pval_assump_LP[1:2,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
boxGood.LP<-apply(Pval_assump_LP[3:4,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
tukGood.LP<-apply(Pval_assump_LP[5:6,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))
logGood.LP<-apply(Pval_assump_LP[7:8,],2,function(x) (x[1]>=0.05)&(x[2]>=0.05))

apply(Pval_assump_LP,1,function(x) sum(is.na(x)))
#bp.meta.LP sw.meta.LP  bp.box.LP  sw.box.LP  bp.tuk.LP  sw.tuk.LP  bp.log.LP  sw.log.LP 
#4           4           4           4           4           4           4           4 

Good.LP=rawGood.LP|logGood.LP|boxGood.LP|tukGood.LP
sum(Good.LP,na.rm = TRUE) #420
sum(!Good.LP,na.rm = TRUE) #228
sum(is.na(Good.LP)) #4
Compliance.LP<-rbind(Good.LP,rawGood.LP, boxGood.LP, tukGood.LP,logGood.LP)
colnames(Compliance.LP)=temp.name

#priority high-low: raw, boxcox, tukey, log
data.choice.LP=vector(length=n.met)
for (i in 1:n.met){
        if(bad.LP[i]){data.choice.LP[i]="noChoice"}
        else if(rawGood.LP[i]==TRUE) {data.choice.LP[i]="raw"}
        else if(boxGood.LP[i]==TRUE) {data.choice.LP[i]="boxcox"}
        else if(tukGood.LP[i]==TRUE) {data.choice.LP[i]="tuk"}
        else if(logGood.LP[i]==TRUE) {data.choice.LP[i]="log"}
        else {data.choice.LP[i]="noChoice"}
}
table(data.choice.LP)
#boxcox      log noChoice      raw      tuk 
# 298        8      232      100       14

bh.LP<-apply(Pval_assump_LP[c(1,3,5,7),],2, function(x) (x[1]>=0.05)|(x[2]>=0.05)|(x[3]>=0.05)|(x[4]>=0.05))
sw.LP<-apply(Pval_assump_LP[c(2,4,6,8),],2, function(x) (x[1]>=0.05)|(x[2]>=0.05)|(x[3]>=0.05)|(x[4]>=0.05))

sum((data.choice.LP=="noChoice")&!sw.LP&!bh.LP,na.rm=TRUE)#17 meet neither, plus 4 bad ==>21
sum((data.choice.LP=="noChoice")&sw.LP&bh.LP,na.rm=TRUE)#24 meet both, but under different transformation

sum((data.choice.LP=="noChoice")&sw.LP,na.rm=TRUE)#127/198 normal but heteroscedastic, use robust errors (HC errors)
sum((data.choice.LP=="noChoice")&bh.LP,na.rm=TRUE)#108/198 homoscedastic but non-normal
sum((data.choice.LP=="noChoice")&bh.LP&!sw.LP,na.rm=TRUE)#84/198 homoscedastic but non-normal, use robust regression
#In conclusion, among the 232 no-choice, 127 with OLS-robust-error, the rest 84 use robust regression, the rest 21 toss.
#OLS-regression for all 420 metabolites that meet both assumptions w/ or w/o transformation.


#Now specify data transformation to use for the 198.The 4 NA metabolites are included but will be properly excluded 
noChoice.HC.LP=(data.choice.LP=="noChoice"&sw.LP) #specify positions for 127 metabolites for HC/robust errors
noChoice.rr.LP=(data.choice.LP=="noChoice"&bh.LP&!sw.LP) #specify positions for the 84 metabolites for robust regression
sum(noChoice.rr.LP&noChoice.HC.LP,na.rm=TRUE) #0, confirm the two don't overlap
#determine which transformation to do, priority high-low: raw, boxcox, tukey, log
#Pval_assump_LP<-rbind(bp.meta.LP,sw.meta.LP,bp.box.LP,sw.box.LP,bp.tuk.LP,sw.tuk.LP,bp.log.LP,sw.log.LP)
data.choice.HC.LP=rep(NA,length=n.met)
for (i in 1:n.met){
        if (!noChoice.HC.LP[i]|bad.LP[i]) next
        if(Pval_assump_LP[2,i]>=0.05) {data.choice.HC.LP[i]="raw"}
        else if(Pval_assump_LP[4,i]>=0.05) {data.choice.HC.LP[i]="boxcox"}
        else if(Pval_assump_LP[6,i]>=0.05) {data.choice.HC.LP[i]="tuk"}
        else if(Pval_assump_LP[8,i]>=0.05) {data.choice.HC.LP[i]="log"}
}
table(data.choice.HC.LP)
#boxcox    log    raw    tuk 
#  69      3     50      5

data.choice.rr.LP=rep(NA,length=n.met)
for (i in 1:n.met){
        if ((!noChoice.rr.LP[i])|bad.LP[i]) next
        if(Pval_assump_LP[1,i]>=0.05) {data.choice.rr.LP[i]="raw"}
        else if(Pval_assump_LP[3,i]>=0.05) {data.choice.rr.LP[i]="boxcox"}
        else if(Pval_assump_LP[5,i]>=0.05) {data.choice.rr.LP[i]="tuk"}
        else if(Pval_assump_LP[7,i]>=0.05) {data.choice.rr.LP[i]="log"}
}
table(data.choice.rr.LP)
#boxcox    log    raw    tuk 
#37      2     43      2 
sum(!is.na(data.choice.rr.LP))#84

#######################################################################
#######(6)CON&LP: regression analysis#################################
#######################################################################

meta.mix.LP=matrix(nrow=n.LP,ncol=n.met)
colnames(meta.mix.LP)=temp.name
dim(meta.mix.LP) #[1]  72 652
meta.mix.LP[,data.choice.LP=="raw"]=meta.LP[,data.choice.LP=="raw"]#cannot use this when vector has na values
meta.mix.LP[,data.choice.LP=="boxcox"]=meta.box.LP[,data.choice.LP=="boxcox"]
meta.mix.LP[,data.choice.LP=="tuk"]=meta.tuk.LP[,data.choice.LP=="tuk"]
meta.mix.LP[,data.choice.LP=="log"]=meta.log.LP[,data.choice.LP=="log"]

meta.mix.LP[,which(data.choice.HC.LP=="raw")]=meta.LP[,which(data.choice.HC.LP=="raw")]#because which() can filter na
meta.mix.LP[,which(data.choice.HC.LP=="boxcox")]=meta.box.LP[,which(data.choice.HC.LP=="boxcox")]
meta.mix.LP[,which(data.choice.HC.LP=="tuk")]=meta.tuk.LP[,which(data.choice.HC.LP=="tuk")]
meta.mix.LP[,which(data.choice.HC.LP=="log")]=meta.log.LP[,which(data.choice.HC.LP=="log")]

meta.mix.LP[,which(data.choice.rr.LP=="raw")]=meta.LP[,which(data.choice.rr.LP=="raw")]
meta.mix.LP[,which(data.choice.rr.LP=="boxcox")]=meta.box.LP[,which(data.choice.rr.LP=="boxcox")]
meta.mix.LP[,which(data.choice.rr.LP=="tuk")]=meta.tuk.LP[,which(data.choice.rr.LP=="tuk")]
meta.mix.LP[,which(data.choice.rr.LP=="log")]=meta.log.LP[,which(data.choice.rr.LP=="log")]

sum(apply(meta.mix.LP,2,function(x) sum(is.na(x))))/n.LP #Confirmed 21 metabolites not added to data for further analysis

#######################################################################
#######(6-a)CON&LP: OLS linear regression for 454 metabolites#########
#######################################################################
library(car)
options(contrasts = c("contr.sum", "contr.poly")) 
test.model=lm(meta.mix.LP[,1]~diet.LP*strain.LP)
summary(test.model)$adj.r.squared #[1] 0.6167251
#http://www.gettinggeneticsdone.com/2011/01/rstats-function-for-extracting-f-test-p.html
#http://blog.minitab.com/blog/adventures-in-statistics-2/what-is-the-f-test-of-overall-significance-in-regression-analysis 
lmp <- function (modelobject) {
        if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
        f <- summary(modelobject)$fstatistic
        p <- pf(f[1],f[2],f[3],lower.tail=F)
        attributes(p) <- NULL
        return(p)
}
lmp(test.model)#[1]  8.480968e-05, this is the pval on the F-test of lm 

test=Anova(test.model, type="III")
test.p<-test[,4]
names(test.p)=row.names(test)# extract pval from ANOVA with names.

OLS.p.LP<-matrix(nrow=length(test.p),ncol=n.met)
colnames(OLS.p.LP)=temp.name
row.names(OLS.p.LP)=row.names(test)

pb <- winProgressBar(title="Anova progress", min=0, max=100, initial=0,label = paste("0% done"))
for(i in 1:n.met){
        if (data.choice.LP[i]=="noChoice") next
        model=lm(meta.mix.LP[,i]~diet.LP*strain.LP)
        test=Anova(model, type="III")
        OLS.p.LP[,i]=test[,4]
        rm(model, test)
        
        info<-sprintf("%d th metabolite analyzed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

#######################################################################
#######(6-b)CON&LP: heteroskedasticity-consistent (HC) errors for 107 metabolites
#######################################################################
#Robust standard errors following lm() can deal with the 107 metabolites in meta.mix.LP[,noChoice.HC.LP] 
#Calculate ANOVA significance and add to the anova.p dataframe.

library(car)
options(contrasts = c("contr.sum", "contr.poly")) 
#http://thestatsgeek.com/2014/02/14/the-robust-sandwich-variance-estimator-for-linear-regression-using-r/
#https://stackoverflow.com/questions/4385436/regression-with-heteroskedasticity-corrected-standard-errors 
#https://stats.stackexchange.com/questions/131401/how-to-get-anova-table-with-robust-standard-errors/132521 ##Wald for robust errors
#http://data.princeton.edu/wws509/notes/c2s3.html ##Basis of Wald test being identical to ANOVA for linear model
#http://data.princeton.edu/wws509/r/robust.html ##Robust standard errors in R
library(sandwich)
library(lmtest)

uracil.LP=lm(meta.mix.LP[,which(temp.name=="uracil")]~diet.LP*strain.LP)
Anova(uracil.LP,type=3)# OLS without HC standard error
data.choice.HC.LP[which(temp.name=="uracil")]#BoxCox

sandwich_se.LP <- diag(vcovHC(uracil.LP, type = "HC1"))^0.5
sandwich_se.LP  ##this is the robust standard errors, same with those in coeftest(, vcov.=vcovHC) ##
coeftest(uracil.LP, vcov. = vcovHC(uracil.LP,type = "HC1"))
Anova(uracil.LP,type=3,vcov. = vcovHC(uracil.LP,type = "HC1"))$`Pr(>F)`
#Coefficient covariances computed by vcovHC(uracil.LP, type = "HC1")
#[1] 7.576170e-01 1.435844e-06 3.652801e-03 5.450389e-01           NA
#[1] "(Intercept)"         "diet.LP"            "strain.LP"          "diet.LP:strain.LP" "Residuals"

row.names(OLS.p.LP) #[1] "(Intercept)"         "diet.LP"            "strain.LP"          "diet.LP:strain.LP" "Residuals"
pb <- winProgressBar(title="Wald with robust error progress", min=0, max=100, initial=0,label = paste("0% done"))
for(i in 1:n.met){
        if ((!noChoice.HC.LP[i])|is.na(noChoice.HC.LP[i])) next
        model=lm(meta.mix.LP[,i]~diet.LP*strain.LP)
        OLS.p.LP[,i]=Anova(model,type=3,vcov. = vcovHC(model,type = "HC1"))$`Pr(>F)`
        rm(model)
        
        info<-sprintf("%d th metabolite analyzed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}
OLS.p.LP[,which(temp.name=="uracil")] ##to confirm loop runned correctly
# (Intercept)            diet.LP          strain.LP diet.LP:strain.LP           Residuals 
# 7.576170e-01      1.435844e-06      3.652801e-03      5.450389e-01                NA 

#######################################################################
#######(6-c)CON&LP: robust regression for 70 metabolites
#######################################################################
library(MASS)
#Mass:rml has issue with list that has unused levels. Convert list to character.
library(car)
options(contrasts = c("contr.sum", "contr.poly")) 

rr.Smethylmethionine <- rlm(meta.mix.LP[,which(temp.name=="S-methylmethionine")]~diet.LP*strain.LP)
summary(rr.Smethylmethionine)
x=linearHypothesis(rr.Smethylmethionine, c("diet.LP1")) 
y=linearHypothesis(rr.Smethylmethionine, c("strain.LP1","strain.LP2","strain.LP3","strain.LP4","strain.LP5","strain.LP6"
                                           ,"strain.LP7")) 
z=linearHypothesis(rr.Smethylmethionine, c("diet.LP1:strain.LP1","diet.LP1:strain.LP2"
                                           ,"diet.LP1:strain.LP3","diet.LP1:strain.LP4"
                                           ,"diet.LP1:strain.LP5","diet.LP1:strain.LP6"
                                           ,"diet.LP1:strain.LP7"))
rr.Smethylmethionine.p=c(x$`Pr(>F)`[2],y$`Pr(>F)`[2],z$`Pr(>F)`[2])
#[1] 2.463804e-95 3.677820e-151 5.283450e-140

test2=lm(meta.mix.LP[,which(temp.name=="S-methylmethionine")]~diet.LP*strain.LP)
Anova(test2,type=3)$`Pr(>F)`        
#[1] 8.406418e-18 5.570605e-01 1.185917e-06 3.009942e-01           NA

#Start a different pvalue table for comparing the two versions for the 70.
sensit.ols.LP=OLS.p.LP
sensit.rr.LP=OLS.p.LP
row.names(sensit.rr.LP)
#[1] "(Intercept)"         "diet.LP"            "strain.LP"          "diet.LP:strain.LP" "Residuals"    

pb <- winProgressBar(title="Sensitivity of the 70 metabolites", min=0, max=100, initial=0,label = paste("0% done"))
for(i in 1:n.met){
        if ((!noChoice.rr.LP[i])|is.na(noChoice.rr.LP[i])) next
        rr.model=rlm(meta.mix.LP[,i]~diet.LP*strain.LP,maxit = 500)
        sensit.rr.LP[2,i]=linearHypothesis(rr.model, c("diet.LP1"))$`Pr(>F)`[2] 
        sensit.rr.LP[3,i]=linearHypothesis(rr.model, c("strain.LP1","strain.LP2","strain.LP3"
                                                       ,"strain.LP4","strain.LP5","strain.LP6"
                                                       ,"strain.LP7"))$`Pr(>F)`[2] 
        sensit.rr.LP[4,i]=linearHypothesis(rr.model, c("diet.LP1:strain.LP1","diet.LP1:strain.LP2"
                                                       ,"diet.LP1:strain.LP3","diet.LP1:strain.LP4"
                                                       ,"diet.LP1:strain.LP5","diet.LP1:strain.LP6"
                                                       ,"diet.LP1:strain.LP7"))$`Pr(>F)`[2]
        ols.model=lm(meta.mix.LP[,i]~diet.LP*strain.LP)
        test=Anova(ols.model, type="III")
        sensit.ols.LP[,i]=test[,4]
        
        rm(rr.model,ols.model, test)
        
        info<-sprintf("%d th metabolite analyzed", i)
        setWinProgressBar(pb, i/n.met*100, label=info)#use %d to print integer, use %% to print a percent sign
}

sum(is.na(sensit.ols.LP[2,which(noChoice.rr.LP)])) #0
sum(is.na(sensit.rr.LP[2,which(noChoice.rr.LP)])) #0
sensit.rr.LP[,which(temp.name=="S-methylmethionine")]
#(Intercept)            diet.LP          strain.LP diet.LP:strain.LP           Residuals 
#       NA      2.463804e-95     3.677820e-151     5.283450e-140                NA 

sensi.diet.LP<-rep(NA,n.met)
sensi.strain.LP<-rep(NA,n.met)
sensi.int.LP<-rep(NA,n.met)
sensi.LP=rbind(sensi.diet.LP,sensi.strain.LP,sensi.strain.LP)
row.names(sensi.LP)=c("diet","strain","diet:strain")
colnames(sensi.LP)=temp.name

for(i in 1:n.met){
        if ((!noChoice.rr.LP[i])|is.na(noChoice.rr.LP[i])) next
        
        if((sensit.ols.LP[2,i]>0.05)&(sensit.rr.LP[2,i]<0.05)){sensi.LP[1,i]="Exclude"} 
        else{sensi.LP[1,i]="Pass"}
        if((sensit.ols.LP[3,i]>0.05)&(sensit.rr.LP[3,i]<0.05)){sensi.LP[2,i]="Exclude"} 
        else{sensi.LP[2,i]="Pass"}
        if((sensit.ols.LP[4,i]>0.05)&(sensit.rr.LP[4,i]<0.05)){sensi.LP[3,i]="Exclude"} 
        else{sensi.LP[3,i]="Pass"}
}

sens.3effects.LP<-apply(sensi.LP,2,function(x) x[1]=="Pass"&x[2]=="Pass"&x[3]=="Pass")
sum(sens.3effects.LP,na.rm=TRUE)
#[1] 33 metabolites passed sensitivity test and the robust regression results for those were kept

#These are metabolites passing sensitivity test with significant diet effect from Robust Regression
pass.sens.diet.rr.LP<-temp.name[which((sens.3effects.LP==TRUE)&sensit.rr.LP[2,]<0.05)]
#[1] "1,2-dipalmitoyl-GPC (16:0/16:0)"             "1-palmitoyl-2-arachidonoyl-GPI (16:0/20:4)*"
#[3] "1-stearoyl-2-linoleoyl-GPA (18:0/18:2)*"     "5-methyluridine (ribothymidine)"            
#[5] "adenosine 5'-diphosphate (ADP)"              "cytidine"                                   
#[7] "gamma-glutamylleucine"                       "kynurenate"                                 
#[9] "N6-carbamoylthreonyladenosine"               "nicotinamide"                               
#[11] "N-linolenoyltaurine *"                       "N-monomethylarginine"                       
#[13] "ophthalmate"                                 "phenylalanylglycine"                        
#[15] "S-adenosylhomocysteine (SAH)"                "taurine"

##Merge p-values from Robust Regression that passes sensitivity test to the master p-value table
Pval.final.LP=OLS.p.LP[2:4,]
row.names(Pval.final.LP)=c("diet.pval","strain.pval", "interaction.pval")
for (i in 1:n.met) {
        if (sens.3effects.LP[i]&!is.na(sens.3effects.LP[i])){
                Pval.final.LP[,i]=sensit.rr.LP[2:4,i]
        }
}

#######################################################################
#######(6-D)CON&LP: Compile regression results########################
#######################################################################
Pval.final.LP[,which(sens.3effects.LP==FALSE)]##Make sure those failed sensitivity all had NA Pval
sum(is.na(Pval.final.LP[1,]))#[1] 72 not analyzed, including 21 not meeting assumptions and 51 failed sensitivity test

Fdr.final.diet.LP=p.adjust(Pval.final.LP[1,],method="fdr", n=sum(!is.na(Pval.final.LP[1,])))#FDR with 592 metabolites
Fdr.final.strain.LP=p.adjust(Pval.final.LP[2,],method="fdr", n=sum(!is.na(Pval.final.LP[2,])))
Fdr.final.int.LP=p.adjust(Pval.final.LP[3,],method="fdr", n=sum(!is.na(Pval.final.LP[3,])))
Fdr.final.LP=rbind(Fdr.final.diet.LP,Fdr.final.strain.LP,Fdr.final.int.LP)
row.names(Fdr.final.LP)=c("diet.fdr","strain.fdr", "interaction.fdr")

colnames(Fdr.final.LP)=temp.name
lm.method.LP=vector(mode="character", length=n.met)
lm.method.LP[which(data.choice.LP!="noChoice")]="OLS-assumptionsMet"
lm.method.LP[which(noChoice.HC.LP)]="HC-error"
lm.method.LP[which(sens.3effects.LP==TRUE)]="RobustRegressionPassingSensitivity"
lm.method.LP[which(sens.3effects.LP==FALSE)]="FailingSensitivity"

table(lm.method.LP)
#                        FailingSensitivity                           HC-error 
# 21                                 51                                127 
#OLS-assumptionsMet RobustRegressionPassingSensitivity 
#420                                 33

lm.transformation.LP=vector(mode="character", length=n.met)
lm.transformation.LP[which(data.choice.LP!="noChoice")]=data.choice.LP[which(data.choice.LP!="noChoice")]
table(lm.transformation.LP)
#     boxcox    log    raw    tuk 
#232    298      8    100     14
lm.transformation.LP[!is.na(data.choice.rr.LP)]=data.choice.rr.LP[!is.na(data.choice.rr.LP)]
table(lm.transformation.LP)
#       boxcox    log    raw    tuk 
# 148    335     10    143     16 
lm.transformation.LP[!is.na(data.choice.HC.LP)]=data.choice.HC.LP[!is.na(data.choice.HC.LP)]
table(lm.transformation.LP)
#       boxcox    log    raw    tuk 
#21    404     13    193     21 



##The compiled results table has 652 columns with each col corresponding to a metabolite##
##Row1-3 have p-values for diet, strain and diet*strain interaction##
##Row4-6 have fdr for diet, strain and diet*strain interaction##
##Row7 tells if data was transformed (raw for no transformation) and which transformation is used (box-cox, tukey or log)##
##Row8 tells which regression method was used:
#### 454 metabolites met both residual normality and residual homoscedasticity, therefore ordinary least squares (OLS) regression and anova to test effects.
#### 107 metabolites only met residual normality, therefore linear regression lwith heteroskedasticity-consistent standard errors and Wald test to test effects.
#### 70 metabolites only met residual homoscedasticity, therefore robust regression and exclusion of 39 metabolites that failed sensitivity test (significance in robust regression but not in OLS)

CompiledResults_LP<-rbind(Pval.final.LP,Fdr.final.LP,lm.transformation.LP,lm.method.LP)
colnames(CompiledResults_LP)=temp.name
#This is the metabolite data table with proper transformation for the regression analysis
Data.wTransformation.LP<-rbind(meta.mix.LP, CompiledResults_LP)
Data.woTransformation.LP<-rbind(meta.LP, CompiledResults_LP)
#Animal sample information is "CC Metabolomics Raw Data File for Analysis_ver62317_ME".

##None of the output files have col names. Use col names of "CC Metabolomics Raw Data File for Analysis_ver62317_ME". 
write.table(Data.wTransformation.LP,"LP_outputFile/LP_regressionResult_dataTransformed.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(Data.woTransformation.LP,"LP_outputFile/LP_regressionResult_dataNotTransformed.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)

##Filter significant diet effects with VIP score
vip=read.delim("OPLSDA&PCA/Output/VIP_foldchange_results.txt", header = F, sep = "\t")
dim(vip)  #[1]   5 653
rowname=vip[2:5,1]
colname=as.character(unlist(vip[1,-1, drop=TRUE]))
vip.new<-sapply(vip[-1,-1],function(x) as.numeric(as.character(x)))
colnames(vip.new)=colname
row.names(vip.new)=rowname
rowname#[1] VIP.new.LVD    foldchange.LVD VIP.new.LP     foldchange.LP
vip.new[1,1]
LVD_regression_VIP=rbind(vip.new[1:2,],Fdr.final.LVD,lm.transformation.LVD,lm.method.LVD)
LP_regression_VIP=rbind(vip.new[3:4,],Fdr.final.LP,lm.transformation.LP,lm.method.LP)
write.table(LVD_regression_VIP,"LVD_outputFile/LVD_regressionResult&VIP.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LP_regression_VIP,"LP_outputFile/LP_regressionResult&VIP.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
LVD_sig.d_VIP.fold1.5=LVD_regression_VIP[,((Fdr.final.diet.LVD<0.05)&!is.na(Fdr.final.diet.LVD))|((vip.new[1,]>1.5)&(abs(vip.new[2,])>1.5))]
LP_sig.d_VIP.fold1.5=LP_regression_VIP[,((Fdr.final.diet.LP<0.05)&!is.na(Fdr.final.diet.LP))|((vip.new[3,]>1.5)&(abs(vip.new[4,])>1.5))]
write.table(LVD_sig.d_VIP.fold1.5,"LVD_outputFile/LVD_regressionResult&VIP_sig&1.5.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LP_sig.d_VIP.fold1.5,"LP_outputFile/LP_regressionResult&VIP_sig&1.5.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)

#############################################################
###OPLSDA analysis for CON-LVD
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
library("ropls")
library(ggplot2)
library(grid)
library(gridExtra)
oplsda_LVD=opls(meta.new.LVD, diet.LVD,predI = 1, orthoI = 1, permI=100)#100 permutations
VIP.LVD<-getVipVn(oplsda_LVD, orthoL = TRUE) #648 metabolites, removed the 4.
LVD.oplscore.p<-getScoreMN(oplsda_LVD, orthoL=FALSE)
LVD.oplscore.o<-getScoreMN(oplsda_LVD, orthoL=TRUE)
LVD.oplscore=cbind(LVD.oplscore.p,LVD.oplscore.o)
LVD.oplscore=as.data.frame(LVD.oplscore)
opl.perc.LVD=c("p1 (4%)", "to1")
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA)
             ,panel.grid.major = element_blank(),panel.grid.minor = element_blank()
             ,strip.background=element_blank(),axis.text.x=element_text(colour="black")
             ,axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black")
             ,plot.margin=unit(c(1,1,1,1),"line"))
p<-ggplot(LVD.oplscore,aes(x=p1,y=o1,color=diet.LVD))
p<-p+geom_point(size=3, alpha=0.5)+theme+xlab(opl.perc.LVD[1]) + ylab(opl.perc.LVD[2]) +
        geom_text(aes(label=strain.LVD),hjust=0, vjust=0) +
        ggtitle("Q2=-0.1, pQ2=0.114")
p
##########################################################################################
#################        Venn diagrams to summarize regression results  ##################
library(gplots)
pdf("venn_lvd_sig.pdf")
op <- par(mar=c(2,2,2,2),mfrow=c(1,1))
sig.d.lvd<-temp.name[(Fdr.final.diet.LVD<0.05)&!is.na((Fdr.final.diet.LVD))]
sig.s.lvd<-temp.name[(Fdr.final.strain.LVD<0.05)&!is.na((Fdr.final.strain.LVD))]
sig.int.lvd<-temp.name[(Fdr.final.int.LVD<0.05)&!is.na((Fdr.final.int.LVD))]
lvd.vip1.5<-temp.name[(vip.new[1,]>1.5)&!is.na((vip.new[1,]))]
lvd.fold1.5<-temp.name[(abs(vip.new[2,])>1.5)&!is.na((vip.new[2,]))]

sig.d.lvd.vip1.5<-intersect(sig.d.lvd,lvd.vip1.5)
sig.d.lvd.fold1.5<-intersect(sig.d.lvd,lvd.fold1.5)

list.lvd<-list(sig.d.lvd, sig.s.lvd, sig.int.lvd)
list.lvd.vipfiltered<-list(sig.d.lvd.vip1.5, sig.s.lvd, sig.int.lvd)
list.lvd.foldfiltered<-list(sig.d.lvd.fold1.5, sig.s.lvd, sig.int.lvd)
list.lvd.d.vip.fold=list(sig.d.lvd,lvd.vip1.5,lvd.fold1.5)

names(list.lvd)=c("d.lvd", "s.lvd", "int.lvd")
names(list.lvd.vipfiltered)=c("d&vip", "s.lvd", "int.lvd")
names(list.lvd.foldfiltered)=c("d&fold", "s.lvd", "int.lvd")
names(list.lvd.d.vip.fold)=c("d.lvd", "vip1.5_lvd","fold1.5_lvd")

venn(list.lvd)
venn(list.lvd.vipfiltered)
venn(list.lvd.foldfiltered)
venn(list.lvd.d.vip.fold)
dev.off()

library(ggplot2)
location.sig.d.lvd=which(Fdr.final.diet.LVD<0.05)
bardata=data.frame(meta.LVD[,location.sig.d.lvd[4]],diet.LVD,strain.LVD)
colnames(bardata)=c("metabolite", "diet", "strain")
p<-ggplot(data=bardata) +
        geom_bar(aes(x=paste(strain,diet,sep="-"), y=metabolite, fill=diet), 
                 position = "dodge", stat = "summary", fun.y = "mean")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        geom_point(data=bardata,aes(x=paste(strain,diet,sep="-"), y=metabolite),
                   position=position_jitter(width =.15), size=0.8) 
p
p1<-ggplot(data=bardata) +
        geom_bar(aes(x=diet, y=metabolite, fill=diet), 
                 position = "dodge", stat = "summary", fun.y = "mean")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        geom_point(data=bardata,aes(x=diet, y=metabolite),
                   position=position_jitter(width =.15), size=1.2) 
p1
p2<-ggplot(data=bardata) +
        geom_bar(aes(x=strain, y=metabolite, fill=strain), 
                 position = "dodge", stat = "summary", fun.y = "mean")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        geom_point(data=bardata,aes(x=strain, y=metabolite),
                   position=position_jitter(width =.15), size=2) 
p2



pdf("venn_lp_sig.pdf")
par(mfrow=c(2,2))
sig.d.lp<-temp.name[(Fdr.final.diet.LP<0.05)&!is.na((Fdr.final.diet.LP))]
sig.s.lp<-temp.name[(Fdr.final.strain.LP<0.05)&!is.na((Fdr.final.strain.LP))]
sig.int.lp<-temp.name[(Fdr.final.int.LP<0.05)&!is.na((Fdr.final.int.LP))]
lp.vip1.5<-temp.name[(vip.new[3,]>1.5)&!is.na((vip.new[3,]))]
lp.fold1.5<-temp.name[(abs(vip.new[4,])>1.5)&!is.na((vip.new[4,]))]

sig.d.lp.vip1.5<-intersect(sig.d.lp,lp.vip1.5)
sig.d.lp.fold1.5<-intersect(sig.d.lp,lp.fold1.5)

list.lp<-list(sig.d.lp, sig.s.lp, sig.int.lp)
list.lp.vipfiltered<-list(sig.d.lp.vip1.5, sig.s.lp, sig.int.lp)
list.lp.foldfiltered<-list(sig.d.lp.fold1.5, sig.s.lp, sig.int.lp)
list.lp.d.vip.fold=list(sig.d.lp,lp.vip1.5,lp.fold1.5)

names(list.lp)=c("d.lp", "s.lp", "int.lp")
names(list.lp.vipfiltered)=c("d&vip", "s.lp", "int.lp")
names(list.lp.foldfiltered)=c("d&fold", "s.lp", "int.lp")
names(list.lp.d.vip.fold)=c("d.lp", "vip1.5_lp","fold1.5_lp")

venn(list.lp)
venn(list.lp.vipfiltered)
venn(list.lp.foldfiltered)
venn(list.lp.d.vip.fold)
dev.off()



c(sum(Fdr.final.diet.LVD<0.05,na.rm = T),sum(Fdr.final.strain.LVD<0.05,na.rm = T),sum(Fdr.final.int.LVD<0.05,na.rm = T))
#[1]   7 496   3
c(sum(Fdr.final.diet.LP<0.05,na.rm = T),sum(Fdr.final.strain.LP<0.05,na.rm = T),sum(Fdr.final.int.LP<0.05,na.rm = T))
#[1] 252 445  79




Meta.pathway<-read.delim("OPLSDA&PCA/InputFile/Metabolite_details_from_raw_metabolonver61617.txt", sep = "\t", header=T)
dim(Meta.pathway)
#[1] 652  13
which(temp.name!=as.character(unlist(Meta.pathway$BIOCHEMICAL)))
#443 [1] "N-acetylglucosaminylasparagine" vs [1] "N-acetylglucosaminylasparagine "
#Have confirmed the metabolites had the same order in all datasets.
LVD.pthw.fdr.vip.fold<-cbind(Meta.pathway, t(LVD_regression_VIP),t(meta.LVD))
LP.pthw.fdr.vip.fold<-cbind(Meta.pathway, t(LP_regression_VIP),t(meta.LP))

LVD.pthw.sigd.vip.fold=LVD.pthw.fdr.vip.fold[(Fdr.final.diet.LVD<0.05)&!is.na(Fdr.final.diet.LVD),]
LP.pthw.sigd.vip.fold=LP.pthw.fdr.vip.fold[(Fdr.final.diet.LP<0.05)&!is.na(Fdr.final.diet.LP),]

write.table(LVD.pthw.sigd.vip.fold,"LVD_outputFile/Pathway_sig.d_VIP&Fold_LVD.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LP.pthw.sigd.vip.fold,"LP_outputFile/Pathway_sig.d_VIP&Fold_LP.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LVD.pthw.fdr.vip.fold,"LVD_outputFile/Pathway_FDR_VIP&Fold_LVD.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LP.pthw.fdr.vip.fold,"LP_outputFile/Pathway_FDR_VIP&Fold_LP.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)

LVD_VIP.fold1.5_fdr=LVD.pthw.fdr.vip.fold[(vip.new[1,]>1.5)&(abs(vip.new[2,])>1.5),]
LP_VIP.fold1.5_fdr=LP.pthw.fdr.vip.fold[(vip.new[3,]>1.5)&(abs(vip.new[4,])>1.5),]
write.table(LVD_VIP.fold1.5_fdr,"LVD_outputFile/LVD_VIP&Fold1.5_wRegression.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(LP_VIP.fold1.5_fdr,"LP_outputFile/LP_VIP&Fold1.5_wRegression.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)

LVD_diet_output=as.character(unlist(diet.LVD))
LVD_diet_output[as.character(unlist(diet.LVD))=="Low Vitamin D"]="LVD"

LP_diet_output=as.character(unlist(diet.LP))
LP_diet_output[as.character(unlist(diet.LP))=="Low Protein"]="LP"

write.csv(rbind(LVD_diet_output,as.character(unlist(strain.LVD))),"LVD_outputFile/LVD_diet.csv",quote=F, row.names = F)
write.csv(rbind(LP_diet_output,as.character(unlist(strain.LP))),"LP_outputFile/LP_diet.csv",quote=F, row.names = F)

All.meta=list(Meta.pathway,LVD.pthw.sigd.vip.fold,LP.pthw.sigd.vip.fold)
enrich.super=lapply(All.meta,function(x) table(x[,3]))
enrich.sub=lapply(All.meta,function(x) table(x[,4]))
names(enrich.super)=c("all","LVD","LP")
names(enrich.sub)=c("all","LVD","LP")

write.table(data.frame(enrich.super),"LVD_outputFile/enrich_super_sig.d.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)
write.table(data.frame(enrich.sub),"LVD_outputFile/enrich_sub_sig.d.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)

Marwa_top15<-read.delim("Marwa_10strainTop15_foldchange_VIP.txt", header=F)
Marwa_top15_vec=rep(NA,15)
Marwa_top15_vec=as.character(Marwa_top15[,1])
Marwa_top15_position=rep(0,15)
for (i in 1:15){
        Marwa_top15_position[i]=which(temp.name==Marwa_top15_vec[i])
}
LVD_vip_fold_regression_MarwaTop15=LVD.pthw.fdr.vip.fold[Marwa_top15_position,]
write.table(LVD_vip_fold_regression_MarwaTop15,"LVD_outputFile/LVD_vip_fold_regression_MarwaTop15.txt", col.names=NA,sep="\t",row.names = TRUE, quote=FALSE)




