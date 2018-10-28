###################################### Variable Selection ###################################
getVarRank = function(los_trdat,crit_indx,finsel_num){
  #importance matrix with each row the importance of one method
  colnum<-length(colnames(los_trdat))-1
  imptc_rec<-rep(0,colnum-2)
  tmpimp_rec<-rep(0,colnum-2)
  tmpdat<-los_trdat[,c(2,3,5:colnum)]
  names(imptc_rec)<-colnames(tmpdat)
  
  trans_dat<-as.data.frame(los_trdat[,2:(colnum+1)])
  trans_dat_scale<-los_trdat[,2:colnum]
  trans_dat_scale$LOS<-log(trans_dat$LOS)
  
  #################### normality check ###############
  fitlm<-lm(LOS~.,data=trans_dat_scale)
  #vif(fitlm)
  #qqnorm(trans_dat_scale$LOS)
  #qqPlot(fitlm,drop.unused.levels = TRUE)

  if(crit_indx==1){
    ############################### importance ranking: model interpretability
  
    ################### robust glm ###########################
    fitrobstlm<-rlm(LOS~.,data=trans_dat_scale,x.ret=TRUE)
    coefval<-fitrobstlm$coefficients
    varimpval<-varImp(fitrobstlm,scale=FALSE)
    ordvarimpval<-varimpval[,1]
    names(ordvarimpval)<-rownames(varimpval)
    ordvarimpval<-ordvarimpval[which(unname(ordvarimpval)>1)]
    ordinfo2<-order(ordvarimpval,decreasing = TRUE)
    rlmvarnam2<-names(ordvarimpval)[ordinfo2]
    if(length(rlmvarnam2)>finsel_num){
      finvarnam<-rlmvarnam2[1:finsel_num]
    }else{
      finvarnam<-rlmvarnam2
    }
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    imptc_rec[updindx]<-1
    
    ####### glm with regulization : Ridge ##################
    col_lis<-unlist(lapply(seq(colnum-1),function(x){ifelse(is.factor(trans_dat_scale[,x]),1,0)}))
    xfactor= model.matrix(LOS~., data=trans_dat_scale[,c(3,which(col_lis==1))])[,-1]
    x<-data.frame(trans_dat_scale[,which(col_lis==0)[-2]],xfactor)
    y= trans_dat_scale$LOS
    names(y)<-'LOS'
    cvfit_ridge<-glmnet(as.matrix(x),as.matrix(y),alpha=0)
    varimpval<-varImp(cvfit_ridge,lambda = 0.1)
    selvarimp<-varimpval[which(varimpval>0),]
    names(selvarimp)<-rownames(varimpval)[which(varimpval>0)]
    selvarind<-order(unname(selvarimp),decreasing = TRUE)
    finvarnam<-names(selvarimp)[selvarind][1:finsel_num]
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### glm with regulization : Lasso ##################
    cvfit_lasso<-glmnet(as.matrix(x),as.matrix(y),alpha=1)
    varimpval<-varImp(cvfit_lasso,lambda = 0.1)
    selvarind<-which(varimpval>0)
    finvarnam<-rownames(varimpval)[selvarind]
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    if(is.vector(continfo)){
      updindx<-which(continfo==TRUE)
    }else{
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    }
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### cox with regulization : Ridge ##################
    col_lis<-unlist(lapply(seq(colnum-1),function(x){ifelse(is.factor(trans_dat[,x]),1,0)}))
    xfactor= model.matrix(LOS~., data=trans_dat[,c(3,which(col_lis==1))])[,-1]
    x<-data.frame(trans_dat[,which(col_lis==0)[-2]],xfactor)
    y= cbind(time=trans_dat$LOS,status=trans_dat$status)
    coxfit_ridge<-glmnet(as.matrix(x), as.matrix(y), family="cox", alpha=0)
    varimpval<-varImp(coxfit_ridge,lambda = 0.1)
    selvarimp<-varimpval[which(varimpval>0),]
    names(selvarimp)<-rownames(varimpval)[which(varimpval>0)]
    selvarind<-order(unname(selvarimp),decreasing = TRUE)
    finvarnam<-names(selvarimp)[selvarind][1:finsel_num]
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### cox with regulization : Lasso ##################
    coxfit_lasso<-glmnet(as.matrix(x), as.matrix(y), family="cox", alpha=1)
    varimpval<-varImp(coxfit_lasso,lambda = 0.1)
    selvarind<-which(varimpval>0)
    finvarnam<-rownames(varimpval)[selvarind]
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### Weibull ##################
    surv_weibfit<-survreg(Surv(LOS,status)~.,data=trans_dat,dist="weibull")
    anvinfo<-anova(surv_weibfit)[-1,]
    anvinfo<-anvinfo[which(anvinfo[,5]<0.05),]
    ordinfo2<-order(anvinfo[,5],decreasing = FALSE)
    varnam<-rownames(anvinfo)[ordinfo2]
    if(length(varnam)>finsel_num){
      finvarnam<-varnam[1:finsel_num]
    }else{
      finvarnam<-varnam
    }
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### Lognormal ##################
    surv_lognfit<-survreg(Surv(LOS,status)~.,data=trans_dat,dist="lognormal")
    anvinfo<-anova(surv_lognfit)[-1,]
    anvinfo<-anvinfo[which(anvinfo[,5]<0.05),]
    ordinfo2<-order(anvinfo[,5],decreasing = FALSE)
    varnam<-rownames(anvinfo)[ordinfo2]
    if(length(varnam)>finsel_num){
      finvarnam<-varnam[1:finsel_num]
    }else{
      finvarnam<-varnam
    }
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### Cox ##################
    coxfit<-coxph(Surv(LOS,status)~.,data=trans_dat)
    anvinfo<-anova(coxfit)[-1,]
    anvinfo<-anvinfo[which(anvinfo[,4]<0.1),]
    ordinfo2<-order(anvinfo[,4],decreasing = FALSE)
    varnam<-rownames(anvinfo)[ordinfo2]
    if(length(varnam)>finsel_num){
      finvarnam<-varnam[1:finsel_num]
    }else{
      finvarnam<-varnam
    }
    continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### relevant feature selection with wrapper algorithm ##
    borsel1<-Boruta(LOS~., data=trans_dat_scale,getImp=getImpRfZ)
    borutavarnam<-getSelectedAttributes(borsel1)
    continfo<-sapply(colnames(tmpdat),grepl,borutavarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    ####### recursive partitioning regression trees #########
    modfit_rp1<-rpart(LOS~., data=trans_dat_scale,method="anova")
    rpartsel<-modfit_rp1$variable.importance
    rpartvarnam<-names(rpartsel)
    continfo<-sapply(colnames(tmpdat),grepl,rpartvarnam)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    modfit_rp2<-rpart(Surv(LOS,status)~.,data=trans_dat,method='exp')
    rpartvarnam2<-names(modfit_rp2$variable.importance)
    continfo<-sapply(colnames(tmpdat),grepl,rpartvarnam2)
    updindx<-apply(continfo,1,function(x){which(x==TRUE)})
    tmpimp_rec[updindx]<-1
    tmpimp_rec[-updindx]<-0
    imptc_rec<-rbind(imptc_rec,tmpimp_rec)
    
    rownames(imptc_rec)<-c('RobustGLM','RidgeGLM','LassoGLM','RidgeCox',
                           'LassoCox','Weibull','Lognormal','Cox',
                           'RFE','RPRT','RPRT_Surv')
  }else{
    
    if(crit_indx==2){
      ############################### stepwise (recursive eliminate): model significance
      ############################## adj R square, AIC etc.
      #glm stepAIC
      fitglm<-glm(LOS~.-1,data=trans_dat_scale)
      glmstep<-stepAIC(fitglm,direction='both',trace=FALSE)
      aicResVarNam<-names(glmstep$coefficients)[-1]
      continfo<-sapply(colnames(tmpdat),grepl,aicResVarNam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      imptc_rec[updindx]<-1
      
      ######################## step Cox
      fitcox<-coxph(Surv(LOS,status)~.,data=trans_dat)
      coxstep<-step(fitcox,direction='both',trace=FALSE)
      aicResVarNam<-names(coxstep$coefficients)
      continfo<-sapply(colnames(tmpdat),grepl,aicResVarNam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      ###########  gradient tree boosting: RMSE #################
      boostguas=gbm(LOS~., data=trans_dat_scale, distribution = "gaussian", interaction.depth = 3, shrinkage = 0.15)
      modsum<-summary(boostguas)[which(summary(boostguas)[,2]>0),]
      if(length(modsum[,1])>finsel_num){
        finvarnam<-as.character(modsum[1:finsel_num,1])
      }else{
        finvarnam<-as.character(modsum[,1])
      }
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      ############ cox RMSE gradient boosting
      boostcox=gbm(Surv(LOS,status)~.,data=trans_dat,distribution='coxph')
      modsum<-summary(boostcox)[which(summary(boostcox)[,2]>0),]
      if(length(modsum[,1])>finsel_num){
        finvarnam<-as.character(modsum[1:finsel_num,1])
      }else{
        finvarnam<-as.character(modsum[,1])
      }
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      ################### decision tree ####################
      tree.model<-tree(LOS~.,data=trans_dat[,-colnum])
      modsum<-summary(tree.model)
      finvarnam<-as.character(modsum$used)
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      X=trans_dat[,-c(3,colnum)]
      y=trans_dat[,3]
      ###### random forest: recursive feature elimination ##
      set.seed(19)
      control<-rfeControl(functions=rfFuncs, method="cv", number=10)
      rfe_modres<- rfe(X, y, sizes=c(8:15), rfeControl=control)
      finvarnam<-predictors(rfe_modres)
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      ######################### random survival forest: minimal depth, importance
      rf_surv<-ranger(Surv(LOS,status) ~ .,data=trans_dat,importance='impurity_corrected')
      varimptval<-rf_surv$variable.importance
      ordindx<-order(abs(unname(varimptval)),decreasing = TRUE)
      selvar_info<-varimptval[ordindx]
      finvarnam<-names(selvar_info)[1:finsel_num]
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      bestsubset_gaus<-bess.one(los_trdat[,-c(1,4)],log(los_trdat[,4]),factor=c(2,4:(colnum-2)),s=10,family='gaussian')
      finvarnam<-names(bestsubset_gaus$bestmodel$coefficients[-1])
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      y<-cbind(trans_dat$LOS,trans_dat$status)
      bestsubset_cox<-bess.one(los_trdat[,-c(1,4)],y,factor=c(2,4:(colnum-2)),s=10,family='cox')
      finvarnam<-names(bestsubset_cox$bestmodel$coefficients)
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      ########### Selection By Filtering ######################
      X=trans_dat[,-c(3,colnum)]
      y=trans_dat[,3]
      rfSBF2 <- rfSBF
      rfSBF2$score <- function(x, y) apply(x, 2, rfSBF$score, y = y)
      sbfCtrl <- sbfControl(functions=rfSBF2,method="repeatedcv",multivariate=TRUE,repeats=5)
      set.seed(14)
      sbfProf <- sbf(X,y,sbfControl = sbfCtrl)
      sbfvarnam<-predictors(sbfProf)
      continfo<-sapply(colnames(tmpdat),grepl,sbfvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
          
      ####################### subset
      bestsub<-regsubsets(LOS~.,data=trans_dat_scale,nvmax=10,nbest=finsel_num)
      modsum<-summary(bestsub)
      finmod<-which.min(modsum$bic)
      finmodmat<-modsum$which[finmod,]
      finvarnam<-names(finmodmat)[which(finmodmat==TRUE)][-1]
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      rownames(imptc_rec)<-c('GLMStep','CoxStep','GradTreeBoost_Gauss','GradTreeBoost_Cox',
                             'DecisionTree','RandomForest_RFE','RandomSurvForest',
                             'BestSubset_Gauss','BestSubset_Cox','SBF','LMSubset')      
    }else{
      ############################# gradient boosting: model discrimination power
      
      ######################### gradient boost for c-index
      boostmodel=gbm(Surv(LOS)~., data=los_trdat[,-1], distribution = list(name='gaussian',metric='conc'),cv.folds=5)
      boostmodel=gbm(Surv(LOS)~., data=los_trdat[,-1], distribution = list(name='coxph',metric='conc'),cv.folds=5)
      finvarnam<-as.character(summary(boostmodel)[1:finsel_num,1])
      continfo<-sapply(colnames(tmpdat),grepl,finvarnam)
      updindx<-apply(continfo,1,function(x){which(x==TRUE)})
      tmpimp_rec[updindx]<-1
      tmpimp_rec[-updindx]<-0
      imptc_rec<-rbind(imptc_rec,tmpimp_rec)
      
      
      boostcox<-glmboost(Surv(LOS)~.,data=los_trdat[,-1],family=CoxPH())
      
      boostcindx<-mboost(Surv(LOS)~.,data=los_trdat[,-1],family=Cindex(sigma=1))
    }
  }

  return(imptc_rec)
}

getCindex_Cox= function(fit_predobj,fitdata){
  cindxval<-rcorr.cens(x=-predict(fit_predobj,newdata=fitdata,type='risk'),S=Surv(fitdata$LOS,fitdata$status))[c(1,3)]
  
  return(cindxval)
}

getCindex_CoxRegu = function(fit_predobj,matx,fitdata){
  predval<-predict.cv.glmnet(fit_predobj,newx=matx,type='response') #link
  cindxval<-rcorr.cens(x=-predval,S=Surv(fitdata$LOS,fitdata$status))[c(1,3)]
  
  return(cindxval)
}

getCindex_Param = function(fitobj,fitdata){
  cindxval<-rcorr.cens(x=predict(fitobj,newdata=fitdata,type='response'),S=Surv(fitdata$LOS,fitdata$status))[c(1,3)]
  #1-survConcordance(Surv(LOS)~predict(fitobj),fitdata)$concordance
  
  return(cindxval)
}

getCindex_GaussianRegu = function(fitobj,matx,fitdata){
  predval<-predict.cv.glmnet(fitobj,newx=matx,s=fitobj$lambda.min,type='response')
  cindxval<-rcorr.cens(x=predval,S=Surv(fitdata$LOS,fitdata$status))[c(1,3)]
  
  return(cindxval)
}

getCindex_CoxTree = function(fitobj,fitdata){
  cindxval<-rcorr.cens(-predict(fitobj,newdata=fitdata,n.trees = 100),Surv(fitdata$LOS,fitdata$status))[c(1,3)]
  
  return(cindxval)
}

getCindex_RF = function(fitobj,fitdata){
  predval<-predict(fitobj,data=fitdata,outcome='test')$predicted.oob
  cindxval<-randomForestSRC::cindex(fitdata$LOS,fitdata$status,predval)
  return(cindxval)
}

getCindex_SVM = function(fitobj,fitdata){
  predval<-predict(fitobj,newdata=fitdata)$predicted
  cindxval<-rcorr.cens(x=predval,S=Surv(fitdata$LOS,fitdata$status))[c(1,3)]
  return(cindxval)
}

#################################### Model Selection #######################################
model_coxreg=function(covdat.train,covdat.test,ret_indx){
  compdata<-as.data.frame(rbind(covdat.train,covdat.test))
  
  if(ret_indx==1){
    modest_res<-coxph(Surv(LOS,status)~.,data=compdata)
    
    return(modest_res)
  }else{
    c_indx<-rep(0,6)
    res.cox_u<-coxph(Surv(LOS,status)~.,data=covdat.train)
    ######## training performance
    c_indx[1:2]<-summary(res.cox_u)$concordance
    ######## testing performance
    c_indx[3:4]<-getCindex_Cox(res.cox_u,covdat.test)
    
    c_indx[5:6]<-getCindex_Cox(res.cox_u,compdata)
    
    return(c_indx)
  }
}

model_regul=function(covdat.train,covdat.test,alphaval,regmod,ret_indx){
  colnum=length(covdat.train[1,])-1
  compdata<-as.data.frame(rbind(covdat.train,covdat.test))
  
  if(ret_indx==1){
    col_lis<-unlist(lapply(seq(2,colnum),function(x){ifelse(is.factor(compdata[,x]),1,0)}))
    xfactor= model.matrix(LOS~., data=compdata[,c(1,which(col_lis==1)+1)])[,-1]
    x<-data.frame(compdata[,which(col_lis==0)+1],xfactor)
    
    set.seed(500+alphaval)
    if(regmod=='cox'){
      y= cbind(time=compdata$LOS,status=compdata$status)
      modest_res<-cv.glmnet(as.matrix(x),as.matrix(y),family=regmod,alpha=alphaval)
    }else{
      y=compdata$LOS
      modest_res<-cv.glmnet(as.matrix(x),as.matrix(y),family=regmod,alpha=alphaval)
    }
    
    return(modest_res)
    
  }else{
    c_indx<-rep(0,6)
    col_lis<-unlist(lapply(seq(2,colnum),function(x){ifelse(is.factor(covdat.train[,x]),1,0)}))
    xfactor= model.matrix(LOS~., data=covdat.train[,c(1,which(col_lis==1)+1)])[,-1]
    x<-data.frame(covdat.train[,which(col_lis==0)+1],xfactor)
    set.seed(100+alphaval)
    col_lisall<-unlist(lapply(seq(2,colnum),function(x){ifelse(is.factor(compdata[,x]),1,0)}))
    if(regmod=='cox'){
      y= cbind(time=covdat.train$LOS,status=covdat.train$status)
      regul_mod<-cv.glmnet(as.matrix(x),as.matrix(y),family=regmod,alpha=alphaval)
      
      ######## training performance
      c_indx[1:2]<-getCindex_CoxRegu(regul_mod,as.matrix(x),covdat.train)
      ######## testing performance
      col_lis<-unlist(lapply(seq(2,colnum),function(x){ifelse(is.factor(covdat.test[,x]),1,0)}))
      xfactor=model.matrix(LOS~., data=covdat.test[,c(1,which(col_lis==1)+1)])[,-1]
      newx<-data.frame(covdat.test[,which(col_lis==0)+1],xfactor)
      c_indx[3:4]<-getCindex_CoxRegu(regul_mod,as.matrix(newx),covdat.test)
      
      xfactor_all=model.matrix(LOS~.,data=compdata[,c(1,which(col_lisall==1)+1)])[,-1]
      newx_all<-data.frame(compdata[,which(col_lisall==0)+1],xfactor_all)
      c_indx[5:6]<-getCindex_CoxRegu(regul_mod,as.matrix(newx_all),compdata)
    }else{
      y=covdat.train$LOS
      regul_mod<-cv.glmnet(as.matrix(x),as.matrix(y),family=regmod,alpha=alphaval)
      ######## training performance
      c_indx[1:2]<-getCindex_GaussianRegu(regul_mod,as.matrix(x),covdat.train)
      ######## testing performance
      xfactor=model.matrix(LOS~., data=covdat.test[,c(1,which(col_lis==1)+1)])[,-1]
      newx<-data.frame(covdat.test[,which(col_lis==0)+1],xfactor)
      c_indx[3:4]<-getCindex_GaussianRegu(regul_mod,as.matrix(newx),covdat.test)
      
      xfactor_all=model.matrix(LOS~.,data=compdata[,c(1,which(col_lisall==1)+1)])[,-1]
      newx_all<-data.frame(compdata[,which(col_lisall==0)+1],xfactor_all)
      c_indx[5:6]<-getCindex_GaussianRegu(regul_mod,as.matrix(newx_all),compdata)
    }
    return(c_indx)
  }
}

model_paramreg = function(covdat.train,covdat.test,diststr,ret_indx){
  compdata=rbind(covdat.train,covdat.test)
  
  if(ret_indx==1){
    modres<-survreg(Surv(LOS,status)~.,data=compdata,dist=diststr)
    
    return(modres)
  }else{
    c_indx<-rep(0,6)
    parmod<-survreg(Surv(LOS,status)~., data=covdat.train, dist=diststr)
    ############# training performance
    c_indx[1:2]=getCindex_Param(parmod,covdat.train)
    ############ testing performance
    c_indx[3:4]=getCindex_Param(parmod,covdat.test)
    c_indx[5:6]=getCindex_Param(parmod,compdata)
    
    return(c_indx)
  }
}

model_dm = function(covdat.train,covdat.test,mod_indx,ret_indx){
  compdata<-rbind(covdat.train,covdat.test)
  
  if(ret_indx==1){
    if(mod_indx==1){
      ############## random forest     survival RF
      modres<-ranger(Surv(LOS,status) ~ .,data=compdata,importance='permutation') #impurity_corrected
    }else{
      if(mod_indx==2){
        ################### boosting tree
        modres<-gbm(LOS~.,data=compdata[,-length(compdata[1,])],distribution='gaussian')
      }else{
        ##################### survival SVM
        modres<-survivalsvm(Surv(LOS,status)~.,data=compdata,type='regression',kernel='poly_kernel',gamma.mu = 0.05)
      }
    } 
    return(modres)    
  }else{
    c_indx<-rep(0,6)
    if(mod_indx==1){      
      rndsurv_forest<-rfsrc(Surv(LOS,status) ~ .,data=covdat.train,mtry = 10, ntree=100)
      c_indx[1]<-getCindex_RF(rndsurv_forest,covdat.train)
      
      ############## testing performance
      c_indx[3]<-getCindex_RF(rndsurv_forest,covdat.test)
      
      c_indx[5]<-getCindex_RF(rndsurv_forest,compdata)
      
    }else{
      if(mod_indx==2){
        ############## boosting tree for COX/GLM?
        modres<-gbm(Surv(LOS,status)~.,data=covdat.train,distribution='coxph')
        ############### training performance
        c_indx[1:2]<-getCindex_CoxTree(modres,covdat.train)
        ############## testing performance
        c_indx[3:4]<-getCindex_CoxTree(modres,covdat.test)
        c_indx[5:6]<-getCindex_CoxTree(modres,compdata)
      }else{
        modres<-survivalsvm(Surv(LOS,status)~.,data=covdat.train,type='regression',kernel="lin_kernel",gamma.mu = 0.1)
        ############### training performance
        c_indx[1:2]<-getCindex_SVM(modres,covdat.train)
        ############## testing performance
        c_indx[3:4]<-getCindex_SVM(modres,covdat.test)
        c_indx[5:6]<-getCindex_SVM(modres,compdata)
      }
    }
    
    return(c_indx)
  }
}


getModCIndex = function(los_trdat,los_tedat){
  cindx_mat<-matrix(0,nrow=1,ncol=6)
  colnames(cindx_mat)<-c('TrainCindex','TrainSe','TestCindex','TestSe','Total','TotalSe')
  
  ############################## Cox regression
  cindx_mat[1,1:6]<-model_coxreg(los_trdat,los_tedat,0)
  ####################### Lasso cox regression
  cindx_mat<-rbind(cindx_mat,model_regul(los_trdat,los_tedat,1,'cox',0))
  ################## Ridge cox regression
  cindx_mat<-rbind(cindx_mat,model_regul(los_trdat,los_tedat,0,'cox',0))
  ####################### Lasso LM regression
  cindx_mat<-rbind(cindx_mat,model_regul(los_trdat,los_tedat,1,'gaussian',0))
  ####################### Ridge regression
  cindx_mat<-rbind(cindx_mat,model_regul(los_trdat,los_tedat,0,'gaussian',0))
  
  ####################### Exponential regression
  cindx_mat<-rbind(cindx_mat,model_paramreg(los_trdat,los_tedat,'exponential',0))  
  ####################### Weibull regression
  cindx_mat<-rbind(cindx_mat,model_paramreg(los_trdat,los_tedat,'weibull',0))
  ####################### Lognormal regression
  cindx_mat<-rbind(cindx_mat,model_paramreg(los_trdat,los_tedat,'lognormal',0))

  ####################### Random forest
  cindx_mat<-rbind(cindx_mat,model_dm(los_trdat,los_tedat,1,0))
  
  ######################## Boosting tree
  cindx_mat<-rbind(cindx_mat,model_dm(los_trdat,los_tedat,2,0))
  ######################### SVM  
  cindx_mat<-rbind(cindx_mat,model_dm(los_trdat,los_tedat,3,0))
  
  rownames(cindx_mat)<-c('Cox','LassoCox','RidgeCox','LassoLM','RidgeLM','Exponential','Weibull',
                         'Lognormal','RF','BoostingTree','SVM')
  
  return(cindx_mat)
}

getModelEst=function(mod_indx,est_trdat,est_tedat){
 modest<-switch(mod_indx,
        model_coxreg(est_trdat,est_tedat,1),
        model_regul(est_trdat,est_tedat,1,'cox',1),
        model_regul(est_trdat,est_tedat,0,'cox',1),
        model_regul(est_trdat,est_tedat,1,'gaussian',1),
        model_regul(est_trdat,est_tedat,0,'gaussian',1),
        model_paramreg(est_trdat,est_tedat,'exponential',1),
        model_paramreg(est_trdat,est_tedat,'weibull',1),
        model_paramreg(est_trdat,est_tedat,'lognormal',1),
        model_dm(est_trdat,est_tedat,1,1),
        model_dm(est_trdat,est_tedat,2,1),
        model_dm(est_trdat,est_tedat,3,1)
 )
 return(modest)
}

############################ generate predictive samples based on estimated model
genPredSamples=function(fitobj,datacomp,mod_type,seedval){
  samplen<-length(datacomp$LOS)
  samp_each=1
  los_predsamp<-rep(0,samplen*samp_each)

  for(i in 1:samplen)
  {
    set.seed(seedval+i)
    if(mod_type==1){
      ################## Cox
      if(datacomp$status[i]==1){
        cur_indvcurv<-survfit(fitobj,newdata = datacomp[i,],censor=FALSE,se.fit = FALSE)
        cur_indvprob = round(runif(samp_each,min=0,max=1),digits=3)
      }else{
        cur_indvcurv<-survfit(fitobj,newdata=datacomp[i,],censor=TRUE,se.fit = FALSE)
        start_prob=cur_indvcurv$surv[which.min(abs(cur_indvcurv$time-datacomp$LOS[i]))]
        cur_indvprob = round(runif(samp_each,min=max(0,start_prob-0.05),max=start_prob+0.01),digits=3)
      }
      cur_los<-c(0,cur_indvcurv$time,max(datacomp$LOS)+1)
      cur_prob<-c(1,cur_indvcurv$surv,0)
            
    }else{
      #Weibull, lognormal and exponential
      cur_los<-predict(fitobj, newdata=datacomp[i,],type="quantile",p=seq(.01,.99,by=.01))
      cur_los<-c(0,cur_los,max(datacomp$LOS)+1)
      cur_prob<-c(1,seq(0.99,0.01,by=-0.01),0)
      if(datacomp$status[i]==1){
        cur_indvprob = round(runif(samp_each,min=0,max=1),digits=3)
      }else{
        start_prob=cur_prob[which.min(abs(cur_los-datacomp$LOS[i]))]
        
        cur_indvprob = round(runif(samp_each,min=max(0,start_prob-0.05),max=start_prob),digits=3)
      }
      
    }
    
    for(k in 1:samp_each){
      cur_diff<-abs(cur_prob-cur_indvprob[k])
      loc_indx<-which.min(cur_diff)
      if(cur_diff[loc_indx]<=1e-4){
        los_predsamp[(i-1)*samp_each+k]=cur_los[loc_indx]
      }else{
        probrange<-NA
        losrange<-NA
        sep_diff<-NA
        if(cur_prob[loc_indx] > cur_indvprob[k]){
          probrange<-c(cur_prob[loc_indx],cur_prob[loc_indx+1])
          losrange<-c(cur_los[loc_indx],cur_los[loc_indx+1])
        }else{
          probrange<-c(cur_prob[loc_indx-1],cur_prob[loc_indx])
          losrange<-c(cur_los[loc_indx-1],cur_los[loc_indx])
        }
        
        if(losrange[1]!=losrange[2]){
          if(probrange[1]!=probrange[2]){
            aprxdat<-approx(losrange,probrange,method='linear',n=round((probrange[1]-probrange[2])*1e4))
            new_loc <- which.min(abs(aprxdat$y-cur_indvprob[k]))
            los_predsamp[(i-1)*samp_each+k]=aprxdat$x[new_loc] #round for int?
          }else{
            los_predsamp[(i-1)*samp_each+k]=(losrange[1]+losrange[2])/2
          }
        }else{
          los_predsamp[(i-1)*samp_each+k]=losrange[1]
        }
      }
    }
  }
  
  return(los_predsamp)
}


######################## Generate samples for boosting tree
genPredSamples_Boost=function(fitobj,datacomp){
  samplen<-length(datacomp$LOS)
  samp_each=1
  los_predsamp<-rep(0,samplen*samp_each)
  
  evalrange=seq(min(datacomp$LOS),max(datacomp$LOS))
  survprob<-matrix(0,nrow=length(datacomp$LOS),length(evalrange))
  for(m in 1:length(evalrange)){
    survprob[,m]<-exp(-basehaz.gbm(datacomp$LOS,
                     datacomp$status,predict(fitobj,newdata=datacomp,n.trees = 100),
                     t.eval=evalrange[m],
                     cumulative = TRUE)*exp(predict(fitobj,newdata=datacomp,n.trees = 100)))
  }
  
  for(i in 1:samplen)
  {
    cur_indvcurv<-survprob[i,]
    cur_los<-c(0,evalrange,365)#max(datacomp$LOS)+1
    cur_prob<-c(1,cur_indvcurv,0)
    
    #inverse sampling
    set.seed(1219+i)
    cur_indvprob = round(runif(samp_each,min=0,max=1),digits=3)
    
    for(k in 1:samp_each){
      cur_diff<-abs(cur_prob-cur_indvprob[k])
      loc_indx<-which.min(cur_diff)
      if(cur_diff[loc_indx]<=1e-4){
        los_predsamp[(i-1)*samp_each+k]=cur_los[loc_indx]
      }else{
        probrange<-NA
        losrange<-NA
        sep_diff<-NA
        if(cur_prob[loc_indx] > cur_indvprob[k]){
          probrange<-c(cur_prob[loc_indx],cur_prob[loc_indx+1])
          losrange<-c(cur_los[loc_indx],cur_los[loc_indx+1])
        }else{
          probrange<-c(cur_prob[loc_indx-1],cur_prob[loc_indx])
          losrange<-c(cur_los[loc_indx-1],cur_los[loc_indx])
        }
        
        if(losrange[1]!=losrange[2]){
          if(probrange[1]!=probrange[2]){
            aprxdat<-approx(losrange,probrange,method='linear',n=round((probrange[1]-probrange[2])*1e4))
            new_loc <- which.min(abs(aprxdat$y-cur_indvprob[k]))
            los_predsamp[(i-1)*samp_each+k]=aprxdat$x[new_loc] #round for int?
          }else{
            los_predsamp[(i-1)*samp_each+k]=(losrange[1]+losrange[2])/2
          }
        }else{
          los_predsamp[(i-1)*samp_each+k]=losrange[1]
        }
      }
    }
  }
  
  return(los_predsamp)
}

######################## Generate samples for survival random forest
genPredSamples_SurvRF=function(fitobj,datacomp){
  samplen<-length(datacomp$LOS)
  samp_each=1
  los_predsamp<-rep(0,samplen*samp_each)
 
  for(i in 1:samplen)
  {
    cur_indvcurv<-predict(fitobj$forest,data=datacomp[i,],outcome='test')
    cur_los<-c(0,cur_indvcurv$unique.death.times,365)#max(datacomp$LOS)+1
    cur_prob<-c(1,cur_indvcurv$survival,0)
    
    #inverse sampling
    set.seed(1219+i)
    cur_indvprob = round(runif(samp_each,min=0,max=1),digits=3)
    
    for(k in 1:samp_each){
      cur_diff<-abs(cur_prob-cur_indvprob[k])
      loc_indx<-which.min(cur_diff)
      if(cur_diff[loc_indx]<=1e-4){
        los_predsamp[(i-1)*samp_each+k]=cur_los[loc_indx]
      }else{
        probrange<-NA
        losrange<-NA
        sep_diff<-NA
        if(cur_prob[loc_indx] > cur_indvprob[k]){
          probrange<-c(cur_prob[loc_indx],cur_prob[loc_indx+1])
          losrange<-c(cur_los[loc_indx],cur_los[loc_indx+1])
        }else{
          probrange<-c(cur_prob[loc_indx-1],cur_prob[loc_indx])
          losrange<-c(cur_los[loc_indx-1],cur_los[loc_indx])
        }
        
        if(losrange[1]!=losrange[2]){
          if(probrange[1]!=probrange[2]){
            aprxdat<-approx(losrange,probrange,method='linear',n=round((probrange[1]-probrange[2])*1e4))
            new_loc <- which.min(abs(aprxdat$y-cur_indvprob[k]))
            los_predsamp[(i-1)*samp_each+k]=aprxdat$x[new_loc] #round for int?
          }else{
            los_predsamp[(i-1)*samp_each+k]=(losrange[1]+losrange[2])/2
          }
        }else{
          los_predsamp[(i-1)*samp_each+k]=losrange[1]
        }
      }
    }
  }
  
  return(los_predsamp)
}


######################### Generate samples for GLMNet or CoxNet
genPredSamples_Regu=function(datacomp){
  samplen<-length(datacomp$LOS)
  samp_each=1
  los_predsamp<-rep(0,samplen*samp_each)
  colnum<-length(datacomp[1,])-1
  col_lis<-unlist(lapply(seq(2,colnum),function(x){ifelse(is.factor(datacomp[,x]),1,0)}))
  xfactor= model.matrix(LOS~., data=datacomp[,c(1,which(col_lis==1)+1)])[,-1]
  covmatx<-data.frame(datacomp[,which(col_lis==0)+1],xfactor)
  datdist=datadist(covmatx)
  y=Surv(datacomp$LOS,datacomp$status)
  options(datadist = 'datdist')
  modest_res = hdcox.lasso(as.matrix(covmatx), y, nfolds = 5, rule = "lambda.min", seed = 1001)
  
  time_dur=seq(min(datacomp$LOS),max(datacomp$LOS))
  surv_curve<-glmnet.survcurve(modest_res$lasso_model,datacomp$LOS, datacomp$status,as.matrix(covmatx),time_dur)$p
  
  for(i in 1:samplen)
  {
    cur_indvcurv<-surv_curve[i,]
    cur_los<-c(0,time_dur,365)
    cur_prob<-c(1,cur_indvcurv,0)
      
    #inverse sampling
    set.seed(1219+i)
    cur_indvprob = round(runif(samp_each,min=0,max=1),digits=3)
    
    for(k in 1:samp_each){
      cur_diff<-abs(cur_prob-cur_indvprob[k])
      loc_indx<-which.min(cur_diff)
      if(cur_diff[loc_indx]<=1e-4){
        los_predsamp[(i-1)*samp_each+k]=cur_los[loc_indx]
      }else{
        probrange<-NA
        losrange<-NA
        sep_diff<-NA
        if(cur_prob[loc_indx] > cur_indvprob[k]){
          probrange<-c(cur_prob[loc_indx],cur_prob[loc_indx+1])
          losrange<-c(cur_los[loc_indx],cur_los[loc_indx+1])
        }else{
          probrange<-c(cur_prob[loc_indx-1],cur_prob[loc_indx])
          losrange<-c(cur_los[loc_indx-1],cur_los[loc_indx])
        }
        
        if(losrange[1]!=losrange[2]){
          if(probrange[1]!=probrange[2]){
            aprxdat<-approx(losrange,probrange,method='linear',n=round((probrange[1]-probrange[2])*1e4))
            new_loc <- which.min(abs(aprxdat$y-cur_indvprob[k]))
            los_predsamp[(i-1)*samp_each+k]=aprxdat$x[new_loc] #round for int?
          }else{
            los_predsamp[(i-1)*samp_each+k]=(losrange[1]+losrange[2])/2
          }
        }else{
          los_predsamp[(i-1)*samp_each+k]=losrange[1]
        }
      }
    }
  }
  
  return(los_predsamp)
}

genPredLOS=function(fitobj,datacomp,seedval,scengrp){
  samplen<-length(datacomp$LOS)
  samp_each=1
  los_predsamp<-rep(0,samplen*samp_each)
  
  for(i in 1:samplen)
  {
    set.seed(seedval+i)
    cur_indvcurv<-survfit(fitobj,newdata = datacomp[i,],censor=FALSE,se.fit = FALSE)
    if(scengrp == 1){
      gengrp = rcat(1,c(0.7,0.15,0.15))
    }else{
      if(scengrp == 2){
        gengrp = rcat(1,c(0.15,0.7,0.15))
      }else{
        gengrp = rcat(1,c(0.15,0.15,0.7))
      }
    }
    
    if(gengrp == 1){
      cur_indvprob = round(runif(samp_each,min=0,max=1),digits=3)
    }else{
      if (gengrp==2){
        cur_indvprob = round(runif(samp_each,min=0,max=0.4),digits=3)
      }else{
        cur_indvprob = round(runif(samp_each,min=0,max=0.2),digits=3)
      }
    }
    
    cur_los<-c(min(datacomp$LOS),cur_indvcurv$time,max(datacomp$LOS)+1)
    cur_prob<-c(1,cur_indvcurv$surv,0)
    
    for(k in 1:samp_each){
      cur_diff<-abs(cur_prob-cur_indvprob[k])
      loc_indx<-which.min(cur_diff)
      if(cur_diff[loc_indx]<=1e-4){
        los_predsamp[(i-1)*samp_each+k]=cur_los[loc_indx]
      }else{
        probrange<-NA
        losrange<-NA
        sep_diff<-NA
        if(cur_prob[loc_indx] > cur_indvprob[k]){
          probrange<-c(cur_prob[loc_indx],cur_prob[loc_indx+1])
          losrange<-c(cur_los[loc_indx],cur_los[loc_indx+1])
        }else{
          probrange<-c(cur_prob[loc_indx-1],cur_prob[loc_indx])
          losrange<-c(cur_los[loc_indx-1],cur_los[loc_indx])
        }
        
        if(losrange[1]!=losrange[2]){
          if(probrange[1]!=probrange[2]){
            aprxdat<-approx(losrange,probrange,method='linear',n=round((probrange[1]-probrange[2])*1e4))
            new_loc <- which.min(abs(aprxdat$y-cur_indvprob[k]))
            los_predsamp[(i-1)*samp_each+k]=aprxdat$x[new_loc] #round for int?
          }else{
            los_predsamp[(i-1)*samp_each+k]=(losrange[1]+losrange[2])/2
          }
        }else{
          los_predsamp[(i-1)*samp_each+k]=losrange[1]
        }
      }
    }
  }
  
  return(los_predsamp)
}
