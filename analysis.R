library(ggplot2)
library(ggthemes)
library(reshape2)
library(readxl)
library(VIM)
library(mice)
library(dplyr)
library(tidyr)
library(psych)
library(ggcorrplot)
library(eRm)
library(ltm)
library(patchwork)
library(difR)
library(semPlot)
library(lavaan)
library(semTools)

#####
#part 1: data preparation, descriptive analyses
#####
{
   
df = read_xlsx("SCS_data.xlsx")
#df = read_xlsx("data.csv")
#if data have been re-downloaded from 
#openpsychometrics, uncomment the above line

SCS_vars = names(df)[1:10]
#set missing values
print(table(df$gender))
df$gender[df$gender == 3] = NA
df[df==0] = NA

print(unique(df$age))
df$age[df$age >= 100] = NA
mean(df$age,na.rm=T)
median(df$age,na.rm=T)
min(df$age,na.rm=T)
max(df$age,na.rm=T)

sprintf("%i cases are incomplete",sum(!complete.cases(df)))
sprintf("%i cases have incomplete SCS data",sum(!complete.cases(df[,SCS_vars])))


#missing data motifs
# and missing proportion per item
pdf("missingplot.pdf",width = 8, height = 4)
aggr(df[!complete.cases(df[,SCS_vars]),SCS_vars], 
     numbers=TRUE, sortVars=TRUE,prop=FALSE,
     labels=SCS_vars, 
     ylab=c("#Cases Missing","Pattern"))
box(which = "figure",lwd=2)
dev.off()

nmissing = rowSums(is.na(df[,SCS_vars]))
table(nmissing[nmissing!=0])
prop.table(table(nmissing[nmissing!=0]))

#missing-at-random analysis
#(check whether missing data points in each variable
#can be jointly predicted by all the other variables)
pvals = data.frame(matrix(ncol = length(SCS_vars), nrow=0))
colnames(pvals) = SCS_vars
for (var in SCS_vars){
  formula = sprintf("I(is.na(%s)) ~ .", var)
  formula0 = sprintf("I(is.na(%s)) ~ 1", var)
  m = summary(glm(formula, data=df[,1:10]))$coefficients
  pvals[var,rownames(m)[2:10]] = m[2:10,"Pr(>|t|)"]
}
min(p.adjust(unlist(pvals), method="fdr"),na.rm=T)


#-> missing at random can be assumed
#remove cases where more than two SCS variables are missing

#15 cases removed
df_clean = df[rowSums(is.na(df[,SCS_vars])) <= 2,]

#use multiple imputation for remaining data
df_clean = complete(mice(df_clean))


#descriptives
df_clean[,1:10] %>% summarise_all(list(mean=mean, median = median,
                                       min = min, max = max)) %>%
  round(1) %>%
  gather(variable, value) %>%
  separate(variable, c("var", "stat"), sep = "\\_") %>%
  spread(var, value) -> descriptives

#fix order of columns in descriptives table
descriptives = descriptives[,c("stat",SCS_vars)]
write.csv(descriptives, "descriptives.csv")
#re-calculate sum score
df_clean$score = rowSums(df_clean[,1:10])

#distribution plot before dichotomization
tmp = melt(
          cbind(data.frame(id=1:nrow(df_clean)),df_clean[,SCS_vars]),
          id.vars="id")
tmp2 = data.frame(table(tmp$variable,tmp$value))
colnames(tmp2) = c("item","response","Freq")
ggplot(tmp2,aes(x=item, y=Freq, fill=response))+geom_col()+theme_clean()
ggsave("distroplot.pdf",width = 4,height = 2)

#dichotomization
dich = df_clean
dich[,1:10] = data.frame(lapply(df_clean[,1:10], 
                                function (x) as.numeric(x > 2)))
dich$score = rowSums(dich[,1:10])

}

#####
#part 2: CTT-style item analysis
#####
{
  
  #tetrachoric correlations
  tetra_cor = tetrachoric(dich[,SCS_vars])
  ggcorrplot(tetra_cor$rho, type = "lower", lab = TRUE)+theme_clean()
  ggsave("tetrachoric_cor_mat.pdf",width = 6, height = 6)
  #dichotomous item statistics (percent and N correct, discriminativity)
  dich.distro = rbind(as.character(round(100*unlist(lapply(dich[,SCS_vars], 
                                                           mean)),1)),
                      as.character(as.integer(unlist(lapply(dich[,SCS_vars], 
                                                            sum)))))
  rownames(dich.distro) = c("item easiness\n(percent in category 1)", 
                            "number of cases in category 1")
  
  discrimination = c()
  for (item in 1:10){
    itemname = SCS_vars[item]
    discrimination[itemname] = as.character(round(biserial(
      rowSums(dich[,-item]),dich[,item]),2))
  }
  
  dich.stats = rbind(dich.distro, discrimination)
  write.csv(dich.stats,"dich_stats.csv")
}


#####
#part 3: estimate and analyze Rasch model
#####

#model fitting
{
  #approach 1: eRm
  #prepare data for eRm estimation
  #(just item data in wide format)
  rasch_model_eRm = RM(dich[,SCS_vars])

  #approach 2: ltm
  #constraint fixes item discriminativity to 1
  rasch_model_ltm = rasch(dich[,SCS_vars],
                          constraint = cbind(length(SCS_vars) + 1, 1))
  smr_ltm = summary(rasch_model_ltm)
  
  

  #TODO check syntax
  #aproach 3: lavaan
  #modified copy from https://jonathantemplin.com/wp-content/uploads/2022/02/
                      #EPSY906_Example05_Binary_IFA-IRT_Models.nb.html
  lavaansyntax = "

    # loadings/discrimination parameters:
    SCS =~ 1*Q1 + 1*Q2 + 1*Q3 + 1*Q4 + 1*Q5 + 1*Q6 + 1*Q7 + 1*Q8 + 1*Q9 + 1*Q10
    
    # threshholds use the | operator and start at value 1 after t:
    Q1 | t1; Q2 | t1; Q3 | t1; Q4 | t1; Q5 | t1; Q6 | t1; Q7 | t1; 
    Q8 | t1; Q9 | t1;Q10 | t1;
    
    # factor mean:
    SCS ~ 0;
    
      # factor variance:
    SCS ~~ 1*SCS
    
    "
  
  rasch_model_lavaan = sem(model = lavaansyntax, data =  dich[,SCS_vars], 
                           ordered = SCS_vars, mimic = "Mplus",
                           estimator = "WLSMV", std.lv = TRUE, 
                           parameterization = "theta")
  smr_lavaan = summary(rasch_model_lavaan, fit.measures = TRUE)
  
  convertTheta2IRT = function(lavObject){
    #modified copy from 
    #https://jonathantemplin.com/wp-content/uploads/2022/02/
        #EPSY906_Example05_Binary_IFA-IRT_Models.nb.html
    
    if (!lavObject@Options$parameterization == "theta") {
      stop("your model is not estimated with parameterization='theta'")
      }
    
    output = inspect(object = lavObject, what = "est")
    if (ncol(output$lambda)>1) { stop("IRT conversion is only valid 
             for one dimensional factor models. 
             Your model has more than one dimension.")
      }    
    a = output$lambda
    b = output$tau/output$lambda
    return(list(a = a, b=b))
  }
  
  #make ICC plot function
  ICC_plot = function(difficulty, discriminativity = 1){
    if (length(discriminativity)==1){
        discriminativity = rep(discriminativity, length(difficulty))
      }
    df = data.frame(x=seq(-6,6,.01))
    for (i in 1:length(difficulty)){
      df[[SCS_vars[i]]] = logistic(x=df$x, d=difficulty[i], 
                                   a=discriminativity[i])
    }
    
    df = melt(df, id.vars = "x")
    colnames(df)[2] = "item"
    plt=ggplot(df, aes(x = x, y = value, color = item, label = item)) + 
      geom_line() + theme_clean() + xlab("Person parameter") + 
      ylab("P(item solved)")
    return(directlabels::direct.label(plt, "last.qp"))
  }
  

  
  #make ICC plots
  difficulties_eRm = -rasch_model_eRm$betapar
  iccplot_eRm=ICC_plot(difficulties_eRm)+ggtitle("eRm")
  

  
  #lme4 difficulties are shifted by .42 from eRm difficulties, why?
  
  difficulties_ltm = smr_ltm$coefficients[1:10,"value"]
  iccplot_ltm = ICC_plot(difficulties_ltm)+ggtitle("ltm")
  
  
  difficulties_lavaan =   convertTheta2IRT(lavObject = rasch_model_lavaan)$b
  
  iccplot_lavaan=ICC_plot(difficulties_lavaan)+ggtitle("lavaan")
  
  
  
  difficulties = rbind( data.frame(model="eRm",
                            item=factor(SCS_vars),
                            difficulty=as.numeric(difficulties_eRm)),
                data.frame(model="ltm",
                           item=factor(SCS_vars),
                           difficulty=as.numeric(difficulties_ltm)),
                data.frame(model="lavaan",
                           item=factor(SCS_vars),
                           difficulty=as.numeric(difficulties_lavaan)),
                data.frame(model="CTT",
                           item=factor(SCS_vars),
                           difficulty=1-as.numeric(dich.distro[1,])/100))
  difficulties_plot = ggplot(difficulties,aes(x=item,y=difficulty,
                                              color=model,group=model)) + 
    geom_point() + geom_line() + theme_clean() + ggtitle("model comparison")+
    scale_x_discrete(breaks=SCS_vars,limits=SCS_vars)
  
  
  difficulties_plot
  ggsave("diffcfig.pdf",width = 4,height = 3)
  
  #arrange plots vertically and save
  iccplot_eRm|iccplot_ltm|iccplot_lavaan

  ggsave("iccfig.pdf",width = 12,height = 3)
  

  #compare fits
    #select second line of output (corresponding to marginal MLE)
  eRm_fit = IC(person.parameter(rasch_model_eRm))[[1]][2,]
  
  ltm_fit = c()
  ltm_fit['value'] = smr_ltm$logLik
  ltm_fit['npar'] = 10
  ltm_fit['AIC'] = smr_ltm$AIC
  ltm_fit['BIC'] = smr_ltm$BIC
  ltm_fit['cAIC'] = NA
  
  rasch_model_fits = rbind(eRm_fit,ltm_fit)
  rownames(rasch_model_fits) = c("eRm","ltm")
  colnames(rasch_model_fits)[1] = "loglik"
  write.csv(rasch_model_fits,"rasch_model_fits.csv")
  
  
  
  #calculate loss per item
  
  predict_responses = function(item_dffc,person_params,item_discr=1){
    item_dffc = as.numeric(item_dffc)
    if (length(item_discr)==1) item_discr = rep(item_discr,length(item_dffc))
    person_params = as.numeric(person_params)
    preds = matrix(nrow=length(person_params),ncol=length(item_dffc))
    for (p in 1:length(person_params)){
      for(i in 1:length(item_dffc)){
        preds[p,i] = logistic(x=person_params[p],d = item_dffc[i], a = item_discr[i])
    }}
    return(preds)
  }
  
  
  #extract latent person abilities
  person_params_eRm = person.parameter(rasch_model_eRm)$theta.table[,"Person Parameter"]
  person_params_ltm=factor.scores(rasch_model_ltm,dich[,SCS_vars])[[1]][,"z1"]
  person_params_lavaan = as.numeric(predict(rasch_model_lavaan))
  
  #make predictions for individual persons and items
  preds_ltm = predict_responses(difficulties_ltm,person_params_ltm)>.5
  preds_eRm = predict_responses(difficulties_eRm,person_params_eRm)>.5
  preds_lavaan = predict_responses(difficulties_lavaan,person_params_lavaan)>.5
  

  #calculate and plot mean 0-1-loss per item
  itemloss_eRm = colMeans(dich[,SCS_vars]!=preds_eRm)
  itemloss_ltm = colMeans(dich[,SCS_vars]!=preds_ltm)
  itemloss_lavaan = colMeans(dich[,SCS_vars]!=preds_lavaan)
  itemloss = rbind(data.frame(model="eRm",item=SCS_vars,loss=itemloss_eRm),
                   data.frame(model="ltm",item=SCS_vars,loss=itemloss_ltm),
                   data.frame(model="lavaan",item=SCS_vars,loss=itemloss_lavaan))
  
  ggplot(itemloss,aes(x=item,y=loss,
                          color=model,group=model)) + 
    geom_point() + geom_line() + theme_clean() + ggtitle("0-1-loss\nmodel comparison")+
    scale_x_discrete(breaks=SCS_vars,limits=SCS_vars)
  
  ggsave("itemlossfig.pdf",width = 4,height = 3)
  
  }

#DIF
{
  data_dif_age = dich[,SCS_vars]
  data_dif_age$age = dich$age > median(dich$age)
  dif_ageL = difLord(data_dif_age,"age",FALSE,"1PL")
  dif_ageR = difRaju(data_dif_age,"age",FALSE, "1PL")
  
  
  data_dif_gender= dich[,c(SCS_vars,"gender")]
  dif_genderL = difLord(data_dif_gender,"gender",1, "1PL")
  dif_genderR = difRaju(data_dif_gender,"gender",1, "1PL")
  
  difstats=data.frame(
    p=-log10(p.adjust(
      c(dif_genderL$p.value,dif_ageL$p.value),method="fdr")),
    item = c(dif_genderL$names,dif_ageL$name),
    groups = c(rep("gender",10),rep("age",10))
    )
  
  
  ggplot(difstats, aes(x=item, y=p, group=groups,col=groups)) + 
    geom_point() + geom_line() + theme_clean() + ylab("-log10(p)")+
    geom_hline(yintercept=-log10(.05))+scale_x_discrete(breaks=SCS_vars,
                                                        limits=SCS_vars)
  
  ggsave("DIF_pvals.pdf",width = 4,height = 3)
  }

#alternative model: 2PL
{
  
  #fit 1PL and 2PL, compare fit
  twoPL_model = ltm(dich[,SCS_vars] ~ z1, IRT.param = TRUE)
  difficulties_2PL = coef(twoPL_model)[,"Dffclt"]
  discriminativities_2PL = coef(twoPL_model)[,"Dscrmn"]
  ICC_2PL = ICC_plot(difficulty = difficulties_2PL,
                     discriminativity = discriminativities_2PL)
  
  
  Rasch_vs_twoPL_comparison = anova(rasch_model_ltm, twoPL_model)
  
  difficulties_1vs2PL = rbind( data.frame(model="Rasch (1-PL)",
                                   item=factor(SCS_vars),
                                   difficulty=as.numeric(difficulties_ltm)),
                        data.frame(model="Birnbaum (2-PL)",
                                   item=factor(SCS_vars),
                                   difficulty=as.numeric(difficulties_2PL)),
                        data.frame(model="CTT",
                                   item=factor(SCS_vars),
                                   difficulty=1-as.numeric(dich.distro[1,])/100))
  ggplot(difficulties_1vs2PL,aes(x=item,y=difficulty,
                                                     color=model,group=model)) + 
    geom_point() + geom_line() + theme_clean() + ggtitle("model comparison")+
    scale_x_discrete(breaks=SCS_vars,limits=SCS_vars)
  ggsave("difficulties_plot_2PL.pdf",width = 6,height = 4)
  
  
  #calculate item-wise infit and outfit
  
  get_outfit = function(ltm_model){
    X=ltm_model$X
    personscores = factor.scores(ltm_model,X)[[1]][,"z1"]
    dffc = coef(ltm_model)[,"Dffclt"]
    discr = coef(ltm_model)[,"Dscrmn"]
    expected = predict_responses(dffc, personscores, discr)
    var_X = expected * (1-expected)
    Z_ij = (X-expected)/sqrt(var_X)
    chisq = colSums(Z_ij**2)
    #divide chisq by n
    return(chisq/nrow(ltm_model$X))
      }
  
  get_infit = function(ltm_model){
    X=ltm_model$X
    personscores = factor.scores(ltm_model,X)[[1]][,"z1"]
    dffc = coef(ltm_model)[,"Dffclt"]
    discr = coef(ltm_model)[,"Dscrmn"]
    expected = predict_responses(dffc, personscores, discr)
    var_X = expected * (1-expected)
    Z_ij = (X-expected)/sqrt(var_X)
    infit = c()
    for (i in 1:length(dffc)){
      infit[i] = sum((var_X[,i] * (Z_ij[,i]**2))/sum(var_X[,i]))}
    return(infit)
  }
  
  outfit_rasch = get_outfit(rasch_model_ltm)
  outfit_2PL = get_outfit(twoPL_model)
  infit_rasch = get_infit(rasch_model_ltm)
  infit_2PL = get_infit(twoPL_model)
  
  inoutfit = rbind(data.frame(model="Rasch",fit="outfit",
                              item=names(outfit_rasch), value=outfit_rasch),
                   data.frame(model="2-PL",fit="outfit",
                              item=names(outfit_2PL), value=outfit_2PL),
                   data.frame(model="Rasch",fit="infit",
                              item=names(outfit_rasch),value=infit_rasch),
                   data.frame(model="2-PL",fit="infit",
                              item=names(outfit_2PL),value=infit_2PL))
  
  ggplot(inoutfit,aes(x=item,y=value,
                                 color=model,group=model)) + facet_wrap(~fit)+
    geom_point() + geom_line() + theme_clean() + ggtitle("model comparison")+
    scale_x_discrete(breaks=SCS_vars,limits=SCS_vars)
  ggsave("inoutfit_plot_2PL.pdf",width = 8,height = 3)
  
}

#polytomous IRT model
{
  grm_constrained = grm(df_clean[,SCS_vars],constrained=T)
  grm_unconstrained = grm(df_clean[,SCS_vars],constrained=F)
  anova(grm_constrained, grm_unconstrained)
  
  #have to save this plot by hand, automatic saving does
  #not work for some reason
  plot(grm_unconstrained, item=2, ask=F)
  
  plot(grm_unconstrained, type="IIC", ask=F,legend=T,cx='topleft',cex=.5)
  
}

#####
#part 4: factorial structure of the data
#####

#factor analyses
{
  
  covdat = cov(df_clean[,SCS_vars])
  N=nrow(df_clean)
  
  
  #unidimensional model
  unidimensional_model <- '
  xi1 =~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  
'
  
  unidimensional_cfa <- cfa(unidimensional_model,
                      sample.cov=covdat,
                      sample.nobs=N,
                      std.lv=T)
  
  
  #correlated traits
  correlated_traits_model <- '
  xi1 =~ Q1+Q2+Q3+Q4+Q10
  xi2 =~ Q5+Q6+Q7+Q8+Q9
  xi1 ~~ xi2
'
  
  correlated_traits_cfa <- cfa(correlated_traits_model,
                               sample.cov=covdat,
                               sample.nobs=N,
                               std.lv=T)
  
  
  #bifactor model (general factor and two item-specific factors)
  bifactor_model <- '
  G =~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  xi1 =~ Q1+Q2+Q3+Q4+Q10
  xi2 =~ Q5+Q6+Q7+Q8+Q9
  G ~~ 0*xi1
  G ~~ 0*xi2
'
  
  bifactor_cfa <- cfa(bifactor_model,
              sample.cov=covdat,
              sample.nobs=N,
              std.lv=T)
  
  
  

  
  #hierarchical model

  hierarchical_model <- '
  xi1 =~ Q1+Q2+Q3+Q4+Q10
  xi2 =~ Q5+Q6+Q7+Q8+Q9
  G =~ xi1+xi2
'
  
  hierarchical_cfa <- cfa(hierarchical_model,
                      sample.cov=covdat,
                      sample.nobs=N,
                      std.lv=T)
  
  # 
  smr_hierarchical = summary(hierarchical_cfa, fit=T)$FIT
  smr_bifactor = summary(bifactor_cfa, fit=T)$FIT
  smr_correlated_traits = summary(correlated_traits_cfa, fit=T)$FIT
  smr_unidimensional = summary(unidimensional_cfa, fit=T)$FIT

  cfaaov_df = data.frame(model=c("hierarchical","bifactor","correlated traits",
                                 "unidimensional"),
                         Df = c(smr_hierarchical["df"],
                                smr_bifactor["df"],
                                smr_correlated_traits["df"],
                                smr_unidimensional["df"]),
                         AIC = c(smr_hierarchical["aic"],
                                 smr_bifactor["aic"],
                                 smr_correlated_traits["aic"],
                                 smr_unidimensional["aic"]),
                         BIC = c(smr_hierarchical["bic"],
                                 smr_bifactor["bic"],
                                 smr_correlated_traits["bic"],
                                 smr_unidimensional["bic"]))

  write.csv(cfaaov_df,"cfaaov_df.csv")



  pdf("semplot_bifactor.pdf", width = 8,height = 4)
  semPaths(bifactor_cfa, "std")
  dev.off()
  
  
  

  
  #fit alternative bifactor model (item Q10 belongs to subscale 2)
  bifactor_model_2 <- '
  G =~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  xi1 =~ Q1+Q2+Q3+Q4
  xi2 =~ Q5+Q6+Q7+Q8+Q9+Q10
  G ~~ 0*xi1
  G ~~ 0*xi2

'
  
  bifactor_cfa_2 <- cfa(bifactor_model_2,
                      sample.cov=covdat,
                      sample.nobs=N,
                      std.lv=T)
  
  #fit alternative bifactor model (item Q10 is its own subscale)
  #-> does not converge
  bifactor_model_3 <- '
  G =~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  xi1 =~ Q1+Q2+Q3+Q4
  xi2 =~ Q5+Q6+Q7+Q8+Q9
  xi3 =~ Q10
  G ~~ 0*xi1
  G ~~ 0*xi2
  G ~~ 0*xi3

'
  
  bifactor_cfa_3 <- cfa(bifactor_model_3,
                        sample.cov=covdat,
                        sample.nobs=N,
                        std.lv=T)
  
  bifactor_model_4 <- '
  G =~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  xi1 =~ Q1+Q2+Q3+Q4
  xi2 =~ Q5+Q6+Q7+Q8+Q9
  xi3 =~ Q10+Q6+Q1
  G ~~ 0*xi1
  G ~~ 0*xi2
  G ~~ 0*xi3

'
  
  bifactor_cfa_4<- cfa(bifactor_model_4,
                        sample.cov=covdat,
                        sample.nobs=N,
                        std.lv=T)
  
  
  pdf("semplot_bifactor_automatic.pdf", width = 8,height = 4)
  semPaths(bifactor_cfa_4, "std")
  dev.off()
  
  
  smr_bif1=summary(bifactor_cfa,fit=T)$FIT
  smr_bif2=summary(bifactor_cfa_2,fit=T)$FIT
  smr_bif3=summary(bifactor_cfa_3,fit=T)$FIT
  smr_bif4=summary(bifactor_cfa_4,fit=T)$FIT
  
  bif_comparison=rbind(smr_bif1,smr_bif2,smr_bif3,smr_bif4)[,c("npar","aic",
                                                               "bic","cfi",
                                                               "tli","srmr",
                                                               "rmsea")]
  rownames(bif_comparison) = c("Bifactor (original)", 
                               "Bifactor (Alternative 1)", 
                               "Bifactor (Alternative 2)",
                               "Bifactor (Automatized)")
  write.csv(bif_comparison, "bifactor_comparison.csv")
  
  }

#reliability, unidimensionality
{
  twofrel = round(semTools::reliability(bifactor_cfa)[-5,],2)
  twofrel = cbind(twofrel,NA)
  twofrel = cbind(twofrel, "two-factor")
  colnames(twofrel) = c("G", "xi_1","xi_2","xi_3","model")
  rownames(twofrel) = c("alpha","omega","omega_2","omega_3")
  
  threefrel = round(semTools::reliability(bifactor_cfa_4)[-5,],2)
  threefrel = cbind(threefrel,"three-factor")
  colnames(threefrel) = colnames(twofrel)
  rownames(threefrel) = rownames(twofrel)
  
  write.csv(rbind(twofrel,threefrel),
            file="composite_reliability.csv")
  

  }

#measurement invariance
{
  
  #fit MI models
  gender_configural = cfa(bifactor_model,
                          data=df_clean[,c(SCS_vars,"gender")],
                          group="gender")
  gender_weak = cfa(bifactor_model,
                          data=df_clean[,c(SCS_vars,"gender")],
                          group="gender",
                          group.equal=c("loadings") )
  gender_strong = cfa(bifactor_model,
                    data=df_clean[,c(SCS_vars,"gender")],
                    group="gender",
                    group.equal=c("loadings","intercepts") )
  gender_strict = cfa(bifactor_model,
                      data=df_clean[,c(SCS_vars,"gender")],
                      group="gender",
                      group.equal=c("loadings","intercepts","residuals") )
  


  mi_gender_modelcomp = anova(gender_configural,
                              gender_weak,
                              gender_strong,
                              gender_strict)
  
  mi_gender_compout=data.frame(cbind(DF=mi_gender_modelcomp$Df,
                          AIC=mi_gender_modelcomp$AIC,
                          BIC=mi_gender_modelcomp$BIC,
        Chisq=mi_gender_modelcomp$Chisq,
        Chisq_diff=mi_gender_modelcomp$`Chisq diff`,
        DF_diff=mi_gender_modelcomp$`Df diff`,
        p=mi_gender_modelcomp$`Pr(>Chisq)`))
  rownames(mi_gender_compout) = c("configural MI",
                                  "weak MI",
                                  "strong MI",
                                  "strict MI")
  write.csv(mi_gender_compout,"mi_gender_compout.csv")
  
  
  
  df_agegroups = df_clean[,c(SCS_vars)]
  df_agegroups$age = df_clean$age > median(df_clean$age)
  age_configural = cfa(bifactor_model_4,
                          data=df_agegroups[,c(SCS_vars,"age")],
                          group="age")
  age_weak = cfa(bifactor_model_4,
                    data=df_agegroups[,c(SCS_vars,"age")],
                    group="age",
                    group.equal=c("loadings") )
  age_strong = cfa(bifactor_model_4,
                      data=df_agegroups[,c(SCS_vars,"age")],
                      group="age",
                      group.equal=c("loadings","intercepts") )
  age_strict = cfa(bifactor_model_4,
                      data=df_agegroups[,c(SCS_vars,"age")],
                      group="age",
                      group.equal=c("loadings","intercepts","residuals") )
  
  mi_age_modelcomp = anova(age_configural,
                              age_weak,
                              age_strong,
                              age_strict)
  
  mi_age_compout=data.frame(cbind(DF=mi_age_modelcomp$Df,
                                     AIC=mi_age_modelcomp$AIC,
                                     BIC=mi_age_modelcomp$BIC,
                                     Chisq=mi_age_modelcomp$Chisq,
                                     Chisq_diff=mi_age_modelcomp$`Chisq diff`,
                                     DF_diff=mi_age_modelcomp$`Df diff`,
                                     p=mi_age_modelcomp$`Pr(>Chisq)`))
  rownames(mi_age_compout) = c("configural MI",
                                  "weak MI",
                                  "strong MI",
                                  "strict MI")
  write.csv(mi_age_compout,"mi_age_compout.csv")
  
  
  #factor structure (just a reminder for
  #modification index interpretation):
  # G =~ Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8+Q9+Q10
  # xi1 =~ Q1+Q2+Q3+Q4
  # xi2 =~ Q5+Q6+Q7+Q8+Q9
  # xi3 =~ Q10+Q6+Q1
  
  #look at mod. indices of strong invariance model
  modindices(gender_weak)
  modindices(age_weak)
  
  
}



