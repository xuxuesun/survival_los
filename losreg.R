library(survival)
library(survminer)
library(mice)
library(car)
library(MASS)
library(leaps)
library(caret)
library(Boruta)
library(rpart)
library(glmnet)
library(gbm)
library(tree)
library(DMwR)
library(gam)
library(ranger)
library(mlbench)
library(BeSS)
library(randomForestSRC)
library(mboost)
library(Hmisc)
library(survAUC)
library(pec)
library(rms)
library(e1071)
library(survivalsvm)
library(c060)
library(hdnom)

library(grid)
library(gridExtra)
library(ggplot2)

###########################################################################################
######################################## Select Covariates Candidates #####################
###########################################################################################

###########################################################################################
######################################## Prepare data #####################################
###########################################################################################
losdat<-read.csv('xxx',stringsAsFactors = FALSE)
therpdat<-read.csv('xxx',stringsAsFactors = FALSE)
covnum=25
losnonna<-losdat[which(complete.cases(losdat$LOS) & (losdat$SpecificAssessmentType != 'Death')),]
covmat<-as.data.frame(matrix(0,nrow=length(losnonna$ptid),ncol=covnum))
colnames(covmat)<-c(colnames(losnonna[,c(1:3,5,7:covnum)]),'StayType','status')

los_rownames<-c(0,rownames(losnonna))

time_invar_orig<-c(1:3,5,8:10)
time_invar_cov<-c(1:4,6:8)
time_var_orig<-c(7,11:covnum)
time_var_cov<-c(5,9:(covnum-2))
for(i in 1:length(losnonna$ptid))
{
  cur_id=losnonna$ptid[i]
  covmat[i,time_invar_cov]<-losnonna[i,time_invar_orig]
  covmat[i,covnum]=1-losnonna$censored[i]
  cur_recind<-which(rownames(losnonna[which(losnonna$ptid==cur_id),])==los_rownames[i+1])
  
  tmpdat<-losdat[seq(as.numeric(los_rownames[i])+1,as.numeric(los_rownames[i+1])),]
  tmpdat<-tmpdat[which(tmpdat$ptid==cur_id),]
  
  if(length(tmpdat$ptid)==1){
    covmat[i,time_var_cov]<-NA
  }else{
    tmpdat2<-tmpdat
    tmpdat2<-tmpdat2[which(tmpdat2$SpecificAssessmentType!='Entry' & tmpdat2$SpecificAssessmentType!='Reentry'),]
    
    if((length(which(tmpdat2$SpecificAssessmentType=='Quarterly'))>0)
       & (covmat$LOS[i]>180)){
      covmat[i,covnum-1]=1   #long term
    }else{
      covmat[i,covnum-1]=0
    }
    
    #initial condition
    if(substring.location(as.character(tmpdat2$SpecificAssessmentType[1]),"DC-RA")$first!=0 |
       substring.location(as.character(tmpdat2$SpecificAssessmentType[1]),"DC-RNA")$first!=0 |
       substring.location(as.character(tmpdat2$SpecificAssessmentType[1]),"DC-EOS")$first!=0)
    {
      covmat[i,time_var_cov]<-NA
      covmat$MedicareStay[i]<-NA
    }else{
      covmat[i,time_var_cov]<-tmpdat2[1,time_var_orig]
    
      curtherapy<-unique(therpdat$therapy[which((therpdat$ptid==tmpdat2[1,]$ptid) &
               (as.character(therpdat$SpecificAssessmentType)==
                  as.character(tmpdat2[1,]$SpecificAssessmentType)))])
      if(length(curtherapy)>1 & cur_recind>1){
        curtherapy<-curtherapy[2]
      }else{
        curtherapy<-curtherapy[1]
      }
      if((covmat$MedicareStay[i]==1) & (curtherapy==1)){
        covmat$MedicareStay[i]=1
      }else{
        covmat$MedicareStay[i]=0
      }
      
    }
  }
}
covmat_cpy<-covmat
covmat<-covmat_cpy[-c(which(is.na(covmat_cpy$TotADL)),which(covmat_cpy$StayType==1)),]

sparsechk<-rep(0,covnum-2)
for(i in 1:(covnum-2)){
  sparsechk[i]<-length(which(covmat[,i]==0))/length(covmat[,i])
}
colnames(covmat)[which(sparsechk>0.9)]

naperc<-rep(0,covnum-2)
for(i in 1:(covnum-2))
{
  naperc[i]<-1-length(na.omit(covmat[,i]))/length(covmat[,i])
}
colnames(covmat)[which(naperc>0.25)]
imputed_Data=mice(as.data.frame(covmat[,which(naperc>0)]), m=5, maxit = 50, method = 'pmm', seed = 500,printFlag=FALSE)
tmp1<-as.matrix(complete(imputed_Data))
covmat[,which(naperc>0)]<-tmp1


comp_regdat<-covmat[which(covmat$LOS<=365),]
comp_regdat<-comp_regdat[,-24]
comp_regdat$CogPerfm_BIMS[which(comp_regdat$CogPerfm_BIMS<8)]<-3 #severe condition
comp_regdat$CogPerfm_BIMS[which(comp_regdat$CogPerfm_BIMS<13 & comp_regdat$CogPerfm_BIMS>=8)]<-2
comp_regdat$CogPerfm_BIMS[which(comp_regdat$CogPerfm_BIMS<=15 & comp_regdat$CogPerfm_BIMS>=13)]<-1

row.names(comp_regdat)<-seq(length(comp_regdat$ptid))

comp_regdat_all<-comp_regdat
comp_regdat<-comp_regdat_all[-remid,]

comp_regdat$Race_Ind<-as.factor(comp_regdat$Race_Ind)
comp_regdat$CogPerfm_BIMS<-as.factor(comp_regdat$CogPerfm_BIMS)
for(i in 6:(length(comp_regdat[1,])-1))
{
  if(length(unique(comp_regdat[,i]))==2){
    comp_regdat[,i]<-as.factor(comp_regdat[,i])
  }
}

comp_regdat<-comp_regdat[,-7]

row.names(comp_regdat)<-seq(length(comp_regdat$ptid))
set.seed(100)
comp_regdat.train = comp_regdat[sample(nrow(comp_regdat),ceiling(0.75*nrow(comp_regdat))),]
comp_regdat.test = comp_regdat[-as.numeric(rownames(comp_regdat.train)),]

source('losreg_funlib.R')
###########################################################################################
######################################## Variable Selection ###############################
###########################################################################################
# importance ranking matrix
finsel_num=10
imp_mat_pval=getVarRank(comp_regdat.train,1,finsel_num)
rank_score_imp<-colSums(imp_mat_pval)
sel_num=length(which(rank_score_imp>=5))#finsel_num
ordindx_imp<-order(unname(rank_score_imp),decreasing = TRUE)
#final selected variabels
finsel_cov_pval<-names(rank_score_imp)[ordindx_imp[1:sel_num]]   

rank_score<-rank_score_imp
ordindx<-ordindx_imp

###########################################################################################
sel_rankvar<-rank_score
names(sel_rankvar)<-c(paste0('X',seq(length(sel_rankvar))))

sel_rankvar<-sel_rankvar[ordindx[1:sel_num]] #score >=5
par(mfrow=c(1,1),mar=c(6.5,5.5,2.5,1.5),family="serif")
barplot(sel_rankvar,main=" ",xlab='Variabels',ylab='Rankings',cex.axis = 1.2,cex.names=1.2,cex.lab=1.5,col=1,names.arg = names(sel_rankvar))

#########################################################################################

continfo<-sapply(colnames(comp_regdat.train),grepl,names(rank_score)[ordindx[1:sel_num]])
updindx<-apply(continfo,1,function(x){which(x==TRUE)})
finsel_covdat.train<-comp_regdat.train[,c(4,updindx,length(comp_regdat.train[1,]))]
finsel_covdat.test<-comp_regdat.test[,c(4,updindx,length(comp_regdat.test[1,]))]
finsel_covdat_all<-comp_regdat[,c(4,updindx,length(comp_regdat[1,]))]

# quick test
coxfit<-coxph(Surv(LOS,status)~.,data=finsel_covdat.train)
summary(coxfit)
logfit<-survreg(Surv(LOS,status)~.,data=finsel_covdat.train,dist='lognormal')
anova(logfit)
coef(logfit)
rcorr.cens(predict(logfit),Surv(finsel_covdat.train$LOS,finsel_covdat.train$status))

###########################################################################################
######################################## Model Selection ##################################
###########################################################################################
cindx_infomat<-getModCIndex(finsel_covdat.train,finsel_covdat.test)
bestmod_indx<-unname(which.max(cindx_infomat[,3]))
bestmod<-names(which.max(cindx_infomat[,3]))

plotres<-data.frame(Model=rownames(cindx_infomat),CindxVal=cindx_infomat[,1],Dataset='Train',stringsAsFactors = FALSE,row.names = NULL)
plotres<-rbind(plotres,data.frame(Model=rownames(cindx_infomat),CindxVal=cindx_infomat[,3],Dataset='Test',stringsAsFactors = FALSE,row.names = NULL))
ggplot(plotres, aes(x=Model, y=CindxVal, fill=Dataset)) + geom_bar(stat = "identity", aes(fill = Dataset), position = position_dodge(width=0.85)) +theme_bw()+ scale_fill_grey(start=0.2,end=0.5) +  labs(x="Model",y="C Index") +ylim(0, 0.9)+ coord_flip()+ geom_hline(yintercept=0.5, linetype="dashed", color = "blue", size = 0.8)+theme(text=element_text(family="serif",size = 11),plot.title = element_text(hjust = 0.5, size = 11),axis.text=element_text(size = 11),legend.text=element_text(size = 11))

plotres1<-data.frame(Model=rownames(cindx_infomat)[1],CindxVal=cindx_infomat[1,1],Dataset='Train')
plotres1<-rbind(plotres1,data.frame(Model=rownames(cindx_infomat)[1],CindxVal=cindx_infomat[1,3],Dataset='Test'))
g1<-ggplot(plotres1, aes(x=Model, y=CindxVal, fill=Dataset)) + geom_bar(stat = "identity", aes(fill = Dataset), width=0.6,position = position_dodge(width=0.5)) +theme_bw()+ scale_fill_grey(start=0.2,end=0.5) +  labs(x="Semi-parametric Hazard",y="C Index") +ylim(0, 0.75)+ coord_flip()+ geom_hline(yintercept=0.6, linetype="dashed", color = "red", size = 0.8)+theme(text=element_text(family="serif",size = 11),plot.title = element_text(hjust = 0.5, size = 11),axis.text=element_text(size = 11),legend.text=element_text(size = 11))

plotres2<-data.frame(Model=rownames(cindx_infomat)[2:5],CindxVal=cindx_infomat[2:5,1],Dataset='Train')
plotres2<-rbind(plotres2,data.frame(Model=rownames(cindx_infomat)[2:5],CindxVal=cindx_infomat[2:5,3],Dataset='Test'))
g2<-ggplot(plotres2, aes(x=Model, y=CindxVal, fill=Dataset)) + geom_bar(stat = "identity", aes(fill = Dataset), position = position_dodge(width=0.85)) +theme_bw()+ scale_fill_grey(start=0.2,end=0.5) +  labs(x="Regularization Based",y="C Index") +ylim(0, 0.75)+ coord_flip()+ geom_hline(yintercept=0.6, linetype="dashed", color = "red", size = 0.8)+theme(text=element_text(family="serif",size = 11),plot.title = element_text(hjust = 0.5, size = 11),axis.text=element_text(size = 11),legend.text=element_text(size = 11))

plotres3<-data.frame(Model=rownames(cindx_infomat)[6:8],CindxVal=cindx_infomat[6:8,1],Dataset='Train')
plotres3<-rbind(plotres3,data.frame(Model=rownames(cindx_infomat)[6:8],CindxVal=cindx_infomat[6:8,3],Dataset='Test'))
g3<-ggplot(plotres3, aes(x=Model, y=CindxVal, fill=Dataset)) + geom_bar(stat = "identity", aes(fill = Dataset), position = position_dodge(width=0.85)) +theme_bw()+ scale_fill_grey(start=0.2,end=0.5) +  labs(x="Parametric Hazard",y="C Index") +ylim(0, 0.75)+ coord_flip()+ geom_hline(yintercept=0.6, linetype="dashed", color = "red", size = 0.8)+theme(text=element_text(family="serif",size = 11),plot.title = element_text(hjust = 0.5, size = 11),axis.text=element_text(size = 11),legend.text=element_text(size = 11))

plotres4<-data.frame(Model=rownames(cindx_infomat)[9:11],CindxVal=cindx_infomat[9:11,1],Dataset='Train')
plotres4<-rbind(plotres4,data.frame(Model=rownames(cindx_infomat)[9:11],CindxVal=cindx_infomat[9:11,3],Dataset='Test'))
g4<-ggplot(plotres4, aes(x=Model, y=CindxVal, fill=Dataset)) + geom_bar(stat = "identity", aes(fill = Dataset), position = position_dodge(width=0.85)) +theme_bw()+ scale_fill_grey(start=0.2,end=0.5) +  labs(x="Data Mining",y="C Index") +ylim(0, 0.75)+ coord_flip()+ geom_hline(yintercept=0.6, linetype="dashed", color = "red", size = 0.8)+theme(text=element_text(family="serif",size = 11),plot.title = element_text(hjust = 0.5, size = 11),axis.text=element_text(size = 11),legend.text=element_text(size = 11))
grid.arrange(g1,g2,g3,g4,nrow=2)

###########################################################################################
######################################## Predicted Samples ################################
###########################################################################################
###################### generate predicted samples for all 4 models
#Cox
cox_est<-getModelEst(1,finsel_covdat.train,finsel_covdat.test)
pval=0
sdval=5
optsd=10
optpval=0
los_coxpred<-genPredSamples(cox_est,finsel_covdat_all,1,sdval)

#Lasso GLM
lassoglm_est<-getModelEst(4,finsel_covdat.train,finsel_covdat.test)
colnum=length(finsel_covdat_all[1,])-1
col_lis<-unlist(lapply(seq(2,colnum),function(x){ifelse(is.factor(finsel_covdat_all[,x]),1,0)}))
xfactor= model.matrix(LOS~., data=finsel_covdat_all[,c(1,which(col_lis==1)+1)])[,-1]
covmatx<-data.frame(finsel_covdat_all[,which(col_lis==0)+1],xfactor)
los_lassogauspred<-as.vector(predict(lassoglm_est,newx=as.matrix(covmatx),type='response',s=lassoglm_est$lambda.min))

excl_col<-c(1,length(finsel_covdat_all[1,]))
#exponential
exp_est<-getModelEst(6,finsel_covdat.train,finsel_covdat.test)
sdval=optsd
los_exppred<-unname(predict(exp_est,newdata=finsel_covdat_all[,-1]))
los_exppred[which(finsel_covdat_all$status==0)]<-genPredSamples(exp_est,finsel_covdat_all[which(finsel_covdat_all$status==0),],2,sdval)

#weibull
weib_est<-getModelEst(7,finsel_covdat.train,finsel_covdat.test)
sdval=optsd
los_weibpred<-unname(predict(weib_est,newdata=finsel_covdat_all[,-1]))
los_weibpred[which(finsel_covdat_all$status==0)]<-genPredSamples(weib_est,finsel_covdat_all[which(finsel_covdat_all$status==0),],2,sdval)

#Lognormal
lognrm_est<-getModelEst(8,finsel_covdat.train,finsel_covdat.test)
sdval=optsd
los_lognrmpred<-unname(predict(lognrm_est,newdata=finsel_covdat_all[,-1]))
los_lognrmpred[which(finsel_covdat_all$status==0)]<-genPredSamples(lognrm_est,finsel_covdat_all[which(finsel_covdat_all$status==0),],2,sdval)

#survival random forest
survrf_est<-getModelEst(9,finsel_covdat.train,finsel_covdat.test)
los_survrf<-genPredSamples_SurvRF(survrf_est,finsel_covdat_all)

# tree GLM
tree_est<-tree(LOS~.,data=finsel_covdat_all[,-length(finsel_covdat_all[1,])])
los_treepred<-unname(predict(tree_est,finsel_covdat_all[,-excl_col]))

#SVM
survsvm_est<-survivalsvm(Surv(LOS,status)~.,data=finsel_covdat_all,type='regression',opt.meth = "quadprog", kernel = "add_kernel",gamma.mu=1)
los_survsvmpred<-as.vector((predict(survsvm_est,newdata=finsel_covdat_all)$predicted))

###########################################################################################
######################################## Performance Evaluation ###########################
###########################################################################################
##################### compare predicted samples and km curve 
##################### for different fitting models (exp, weibull, ml)
survPlots=list()
survobj_ref<-Surv(finsel_covdat_all$LOS,finsel_covdat_all$status)
reffit<-survfit(survobj_ref~1)

survobj_cox<-Surv(los_coxpred,finsel_covdat_all$status)
modfit_cox<-survfit(survobj_cox~1)
survPlots[[1]]<-ggsurvplot_combine(list(reffit,modfit_cox),
            data=finsel_covdat_all, title='Cox',xlab = "Length of Stay",censor=FALSE,
            ylab='Stay-in-NH Probability',legend.title='',
            palette=c('red','blue'),
            legend.labs=c("Real Data","Predicted Samples"),
            conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
              theme_update(panel.grid.minor = element_line(colour='grey'),
                legend.key=element_blank(),
                           legend.background=element_blank(),
                           legend.box.background=element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                           plot.title = element_text(hjust = 0.5)))

survobj_lassogaus<-Surv(los_lassogauspred)
modfit_lassogaus<-survfit(survobj_lassogaus~1)
survPlots[[2]]<-ggsurvplot_combine(list(reffit,modfit_lassogaus),
          data=finsel_covdat_all,title='Lasso GLM',xlab="Length of Stay",censor=FALSE,
          ylab='Stay-in-NH Probability',legend.title='',
          palette=c('red','blue'),
          legend.labs=c("Real Data","Predicted Samples"),
          conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
            theme_update(panel.grid.minor = element_line(colour='grey'),
                         legend.key=element_blank(),
                         legend.background=element_blank(),
                         legend.box.background=element_blank(),
                         panel.background = element_blank(), axis.line = element_line(colour = "black"),
                         plot.title = element_text(hjust = 0.5)))

survobj_exp<-Surv(los_exppred,finsel_covdat_all$status)
modfit_exp<-survfit(survobj_exp~1)
survPlots[[3]]<-ggsurvplot_combine(list(reffit,modfit_exp),
    data=finsel_covdat_all,title='Exponential', xlab = "Length of Stay",censor=FALSE,
    ylab='Stay-in-NH Probability',legend.title='',
    palette=c('red','blue'),
    legend.labs=c("Real Data","Predicted Samples"),
    conf.int = FALSE,override.aes=list(shape=NA),ggtheme = theme_pubclean(base_family = 'Times New Roman')+
      theme_update(panel.grid.minor = element_line(colour='grey'),
                   legend.key=element_blank(),
                   legend.background=element_blank(),
                   legend.box.background=element_blank(),
                   panel.background = element_blank(), axis.line = element_line(colour = "black"),
                   plot.title = element_text(hjust = 0.5)))

survobj_weib<-Surv(los_weibpred,finsel_covdat_all$status)
modfit_weib<-survfit(survobj_weib~1)
survPlots[[4]]<-ggsurvplot_combine(list(reffit,modfit_weib),
  data=finsel_covdat_all,title='Weibull', xlab = "Length of Stay",censor=FALSE,
  ylab='Stay-in-NH Probability',legend.title='',
  palette=c('red','blue'),
  legend.labs=c("Real Data","Predicted Samples"),
  conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
    theme_update(panel.grid.minor = element_line(colour='grey'),
                 legend.key=element_blank(),
                 legend.background=element_blank(),
                 legend.box.background=element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 plot.title = element_text(hjust = 0.5)))

survobj_lognrm<-Surv(los_lognrmpred,finsel_covdat_all$status)
modfit_lognrm<-survfit(survobj_lognrm~1)
survPlots[[5]]<-ggsurvplot_combine(list(reffit,modfit_lognrm),
  data=finsel_covdat_all,title='Lognormal', xlab = "Length of Stay",censor=FALSE,
  ylab='Stay-in-NH Probability',legend.title='',
  palette=c('red','blue'),
  legend.labs=c("Real Data","Predicted Samples"),
  conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
    theme_update(panel.grid.minor = element_line(colour='grey'),
                 legend.key=element_blank(),
                 legend.background=element_blank(),
                 legend.box.background=element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 plot.title = element_text(hjust = 0.5)))

survobj_survrf<-Surv(los_survrf,finsel_covdat_all$status)
modfit_survrf<-survfit(survobj_survrf~1)
survPlots[[6]]<-ggsurvplot_combine(list(reffit,modfit_survrf),
      data=finsel_covdat_all,title='Random Survival Forest',censor=FALSE,
      xlab = "Length of Stay",ylab='Stay-in-NH Probability',legend.title='',
      palette=c('red','blue'),
      legend.labs=c("Real Data","Predicted Samples"),
      conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
        theme_update(panel.grid.minor = element_line(colour='grey'),
                     legend.key=element_blank(),
                     legend.background=element_blank(),
                     legend.box.background=element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"),
                     plot.title = element_text(hjust = 0.5)))

survobj_tree<-Surv(los_treepred)
modfit_tree<-survfit(survobj_tree~1)
survPlots[[7]]<-ggsurvplot_combine(list(reffit,modfit_tree),
  data=finsel_covdat_all,title='Decision Tree',
  xlab = "Length of Stay",ylab='Stay-in-NH Probability',
  legend.title='',censor=FALSE,
  palette=c('red','blue'),
  legend.labs=c("Real Data","Predicted Samples"),
  conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
    theme_update(panel.grid.minor = element_line(colour='grey'),
                 legend.key=element_blank(),
                 legend.background=element_blank(),
                 legend.box.background=element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"),
                 plot.title = element_text(hjust = 0.5)))

survobj_svm<-Surv(los_survsvmpred,finsel_covdat_all$status)
modfit_svm<-survfit(survobj_svm~1)
survPlots[[8]]<-ggsurvplot_combine(list(reffit,modfit_svm),
   data=finsel_covdat_all,title='Survival SVM',
   xlab = "Length of Stay",ylab='Stay-in-NH Probability',
   legend.title='',censor=FALSE,
   palette=c('red','blue'),
   legend.labs=c("Real Data","Predicted Samples"),
   conf.int = FALSE,ggtheme = theme_pubclean(base_family = 'Times New Roman')+
     theme_update(panel.grid.minor = element_line(colour='grey'),
                  legend.key=element_blank(),
                  legend.background=element_blank(),
                  legend.box.background=element_blank(),
                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                  plot.title = element_text(hjust = 0.5)))

grid.arrange(survPlots[[1]]$plot,survPlots[[2]]$plot,survPlots[[3]]$plot,
             survPlots[[4]]$plot,survPlots[[5]]$plot,survPlots[[6]]$plot,
             survPlots[[7]]$plot,survPlots[[8]]$plot,ncol = 4, nrow = 2)   #,layout_matrix=laymat


ksvalmat<-matrix(0,nrow=8,ncol=2)
colnames(ksvalmat)<-c('D','p value')
rownames(ksvalmat)<-c('Cox','Lasso','Exponential',
                      'Weibull','Lognormal','RF','DT','SVM')

los_predmat<-rbind(los_coxpred,los_lassogauspred,los_exppred,
                   los_weibpred,los_lognrmpred,los_survrf,los_treepred,los_survsvmpred)

for(i in 1:8){
  tmpks=ks.test(finsel_covdat_all$LOS,los_predmat[i,]) #jitter to bypass ties issue
  ksvalmat[i,]=c(unname(tmpks$statistic),tmpks$p.value)
}

logrkvalmat<-matrix(0,ncol=2,nrow=8)
colnames(logrkvalmat)<-c('Chisq','pvalue')
rownames(logrkvalmat)<-rownames(ksvalmat)
for(i in 1:8){
  if(i == 2 | i==7){
    cur_dat<-as.data.frame(rbind(cbind(finsel_covdat_all$LOS,finsel_covdat_all$status,label=rep(0,length(finsel_covdat_all$LOS))),
     cbind(los_predmat[i,],rep(1,length(los_predmat[i,])),label=rep(1,length(los_predmat[i,])))))
    rhoval=1
  }else{
    cur_dat<-as.data.frame(rbind(cbind(finsel_covdat_all$LOS,finsel_covdat_all$status,label=rep(0,length(finsel_covdat_all$LOS))),
                                 cbind(los_predmat[i,],finsel_covdat_all$status,label=rep(1,length(los_predmat[i,])))))
    rhoval=0
  }
  colnames(cur_dat)<-c('LOS','status','label')
  
  survdiffval<-survdiff(formula = Surv(LOS, status) ~ label, data = cur_dat,rho=rhoval)
  
  logrkvalmat[i,1]<-survdiffval$chisq
  logrkvalmat[i,2]<-pchisq(survdiffval$chisq, length(survdiffval$n)-1, lower.tail = FALSE)
}


############################################################
############################################################
#################### link trial #########################
########################### LoS generation ##############
covset_dt<-data.table(finsel_covdat_all[,-c(1,6,7,11)])
dup_rec<-duplicated(covset_dt,by=colnames(covset_dt))
uniq_covset1=unique(covset_dt[dup_rec])
totuniq<-rbind(unique(covset_dt,by=colnames(covset_dt)),uniq_covset1)
scen_cov1<-matrix(0,nrow=uniq_covset,ncol=finsel_num-2)


################## covariates effect ######################
########### group 1
library(extraDistr)
los_grp1<-rep(0,1000)
grp_dat1<-as.data.frame(matrix(0,nrow=1000,ncol=7))
grp_dat1[,1]<-c(rep(max(finsel_covdat_all$LOS),500),rep(min(finsel_covdat_all$LOS),500))
colnames(grp_dat1)<-colnames(finsel_covdat_all)[-c(2,3,6,7)]
set.seed(30)
adlgrp1 = rcat(1000,c(0.7,0.15,0.15))
for(i in 1:1000){
  grp_dat1[i,7]=1
  set.seed(300+(i-1)*5+1)
  grp_dat1[i,4]=rcat(1,c(0.15,0.85))-1
  
  set.seed(300+(i-1)*5+3)
  if(adlgrp1[i]==1){
    grp_dat1[i,2]=rdunif(1,0,4)
  }else{
    if(adlgrp1[i]==2){
      grp_dat1[i,2]=rdunif(1,5,10)
    }else{
      grp_dat1[i,2]=rdunif(1,11,16)
    }
  }
  
  set.seed(300+(i-1)*5+2)
  if(adlgrp1[i]==1){
    grp_dat1[i,3]=rdunif(1,39,55)
  }else{
    if(adlgrp1[i]==2){
      grp_dat1[i,3]=rdunif(1,56,90)
    }else{
      grp_dat1[i,3]=rdunif(1,90,102)
    }
  }
    
  set.seed(300+(i-1)*5+4)
  grp_dat1[i,5]=rcat(1,c(0.8,0.1,0.1))
  set.seed(300+(i-1)*5+5)
  grp_dat1[i,6]=rcat(1,c(0.85,0.15))-1
}

par(mfrow=c(2,2),mar=c(6.5,5.5,2.5,1.5),family="serif")
hist(grp_dat1[,2],breaks=16,xlab='Total ADL',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat1[,3],xlab='Age',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat1[,5],xlab='Cognitive Performance',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat1[,6],xlab='Psychological disease',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)


grp_dat1[,4]<-as.factor(grp_dat1[,4])
grp_dat1[,5]<-as.factor(grp_dat1[,5])
grp_dat1[,6]<-as.factor(grp_dat1[,6])

############## check individual curve
los_grp1 = genPredLOS(tmp1,grp_dat1,705,1)
los_grp1 = round(los_grp1)
par(mfrow=c(1,1),mar=c(6.5,5.5,2.5,1.5),family="serif")
hist(los_grp1,xlab='LOS',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)

set.seed(80)
adlgrp2 = rcat(1000,c(0.25,0.5,0.25))
los_grp2<-rep(0,1000)
grp_dat2<-as.data.frame(matrix(NA,nrow=1000,ncol=7))
grp_dat2[,1]<-c(rep(max(finsel_covdat_all$LOS),500),rep(min(finsel_covdat_all$LOS),500))
colnames(grp_dat2)<-colnames(finsel_covdat_all)[-c(2,3,6,7)]
for(i in 1:1000){
  grp_dat2[i,7]=1
  set.seed(2500+(i-1)*5+1)
  grp_dat2[i,4]=rcat(1,c(0.5,0.5))-1
  
  set.seed(2500+(i-1)*5+2)
  if(adlgrp2[i]==1){
    grp_dat2[i,3]=rdunif(1,39,55)
  }else{
    if(adlgrp2[i]==2){
      grp_dat2[i,3]=rdunif(1,56,90)
    }else{
      grp_dat2[i,3]=rdunif(1,90,102)
    }
  }
  
  set.seed(2500+(i-1)*5+3)
  if(adlgrp2[i]==1){
    grp_dat2[i,2]=rdunif(1,0,4)
  }else{
    if(adlgrp2[i]==2){
      grp_dat2[i,2]=rdunif(1,5,10)
    }else{
      grp_dat2[i,2]=rdunif(1,11,16)
    }
  }
  
  set.seed(2500+(i-1)*5+4)
  grp_dat2[i,5]=rcat(1,c(0.25,0.5,0.25))
  set.seed(2500+(i-1)*5+5)
  grp_dat2[i,6]=rcat(1,c(0.5,0.5))-1
}

par(mfrow=c(2,2),mar=c(6.5,5.5,2.5,1.5),family="serif")
hist(grp_dat2[,2],breaks=16,xlab='Total ADL',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat2[,3],xlab='Age',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat2[,5],xlab='Cognitive Performance',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat2[,6],xlab='Psychological disease',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)

############# check individual curve determine conditional prob

grp_dat2[,4]<-as.factor(grp_dat2[,4])
grp_dat2[,5]<-as.factor(grp_dat2[,5])
grp_dat2[,6]<-as.factor(grp_dat2[,6])

los_grp2 = genPredLOS(tmp1,grp_dat2,705,2)
los_grp2 = round(los_grp2)
par(mfrow=c(1,1),mar=c(6.5,5.5,2.5,1.5),family="serif")
hist(los_grp2,xlab='LOS',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)


set.seed(120)
adlgrp3 = rcat(1000,c(0.15,0.15,0.7))
los_grp3<-rep(0,1000)
grp_dat3<-as.data.frame(matrix(NA,nrow=1000,ncol=7))
grp_dat3[,1]<-c(rep(max(finsel_covdat_all$LOS),500),rep(min(finsel_covdat_all$LOS),500))
colnames(grp_dat3)<-colnames(finsel_covdat_all)[-c(2,3,6,7)]
for(i in 1:1000){
  grp_dat3[i,7]=1
  set.seed(5000+(i-1)*5+1)
  grp_dat3[i,4]=rcat(1,c(0.85,0.15))-1
  
  set.seed(5000+(i-1)*5+2)
  if(adlgrp3[i]==1){
    grp_dat3[i,3]=rdunif(1,39,55)
  }else{
    if(adlgrp3[i]==2){
      grp_dat3[i,3]=rdunif(1,56,90)
    }else{
      grp_dat3[i,3]=rdunif(1,90,102)
    }
  }
  
  set.seed(5000+(i-1)*5+3)
  if(adlgrp3[i]==1){
    grp_dat3[i,2]=rdunif(1,0,4)
  }else{
    if(adlgrp3[i]==2){
      grp_dat3[i,2]=rdunif(1,5,10)
    }else{
      grp_dat3[i,2]=rdunif(1,11,16)
    }
  }
  
  set.seed(5000+(i-1)*5+4)
  grp_dat3[i,5]=rcat(1,c(0.1,0.1,0.8))
  set.seed(2500+(i-1)*5+5)
  grp_dat3[i,6]=rcat(1,c(0.15,0.85))-1
}

par(mfrow=c(2,2),mar=c(6.5,5.5,2.5,1.5),family="serif")
hist(grp_dat3[,2],breaks=16,xlab='Total ADL',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat3[,3],xlab='Age',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat3[,5],xlab='Cognitive Performance',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)
hist(grp_dat3[,6],xlab='Psychological disease',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)

############# check individual curve determine conditional prob
grp_dat3[,4]<-as.factor(grp_dat3[,4])
grp_dat3[,5]<-as.factor(grp_dat3[,5])
grp_dat3[,6]<-as.factor(grp_dat3[,6])

los_grp3 = genPredLOS(tmp1,grp_dat3,705,3)
los_grp3 = round(los_grp3)
par(mfrow=c(1,1),mar=c(6.5,5.5,2.5,1.5),family="serif")
hist(los_grp3,xlab='LOS',ylab='Frequency',main=' ',col='blue',cex.lab=2,cex.axis=2)

hist(los_grp1, xlim=c(0,200), col=rgb(1,0,1,0.5), xlab="LOS",
     ylab='Frequency', main=" ",cex.lab=2,cex.axis=2)
hist(los_grp2,xlim=c(0,200), col=rgb(1,1,0,0.7), add=T,cex.lab=2,cex.axis=2)
hist(los_grp3,xlim=c(0,200),col=rgb(0,1,0,0.5),add=T,cex.lab=2,cex.axis=2)
legend("topright", legend=c("Group1","Group2",'Group3'), col=c(rgb(1,0,1,0.5),
     rgb(1,1,0,0.7),rgb(0,1,0,0.5)),pt.cex=2, pch=15)


lostab<-rbind(los_grp1,los_grp2,los_grp3)
write.table(lostab,file='xxx',sep=',',row.names = FALSE,col.names = FALSE)