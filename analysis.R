require(ggplot2)
require(ggthemes)
require(reshape2)
require(readxl)
require(VIM)
require(mice)
require(dplyr)
require(tidyr)
require(psych)
require(ggcorrplot)
require(eRm)
require(ltm)
require(lavaan)
require(patchwork)
require(difR)

#####
#part 1: data preparation, descriptive analyses
#####
{
df = read_xlsx("SCS_data.xlsx")
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
  
  #biserial correlations
  biserial_cor = biserial(dich[,SCS_vars],dich[,SCS_vars])
  ggcorrplot(biserial_cor, type = "lower", lab = TRUE)+theme_clean()
  ggsave("biserial_cor_mat.pdf",width = 6, height = 6)
  #dichotomous item statistics (percent and N correct, discriminativity)
  dich.distro = rbind(as.character(round(100*unlist(lapply(dich[,SCS_vars], 
                                                           mean)),1)),
                      as.character(as.integer(unlist(lapply(dich[,SCS_vars], 
                                                            sum)))))
  rownames(dich.distro) = c("percent in category 1", 
                            "number of cases in category 1")
  
  discrimination = c()
  for (item in 1:10){
    itemname = SCS_vars[item]
    discrimination[itemname] = as.character(round(biserial(
      rowSums(dich[,-item]),dich[,item]),2))
  }
  
  dich.stats = rbind(dich.distro, discrimination)
}


#####
#part 2: estimate Rasch model
#####
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
  smr_lavaan = summary(rasch_model_lavaan, fit.measures = TRUE,
                       rsquare = TRUE, standardized = TRUE)
  
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
  
  #TODO: check ICC plotting fct
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
    scale_x_discrete(breaks=paste0("Q",1:10),limits=paste0("Q",1:10))
  
  
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
  
  fitdf = rbind(eRm_fit,ltm_fit)
  rownames(fitdf) = c("eRm","ltm")
  
  #TODO: look up good ranges of lavaan models
  smr_lavaan$FIT
  
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
    geom_hline(yintercept=-log10(.05))+scale_x_discrete(breaks=paste0("Q",1:10),
                                                        limits=paste0("Q",1:10))
  
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
                                   difficulty=as.numeric(dich.distro[1,])/100))
  ggplot(difficulties_1vs2PL,aes(x=item,y=difficulty,
                                                     color=model,group=model)) + 
    geom_point() + geom_line() + theme_clean() + ggtitle("model comparison")
  ggsave("difficulties_plot_2PL.pdf",width = 6,height = 8)
  
  }

#alternative models: bifactor, ...
{
  
}

#reliability, unidimensionality
{
  
  
}

#polytomous IRT model
{
  grm_model = grm(df_clean[,SCS_vars],constrained=T)
  plot(grm_model)
}