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
require(lme4)
require(lavaan)


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
df_clean[,1:10] %>% summarise_all(list(mean=mean, median = median, min = min, max = max)) %>%
  round(1) %>%
  gather(variable, value) %>%
  separate(variable, c("var", "stat"), sep = "\\_") %>%
  spread(var, value) -> descriptives

#re-calculate sum score
df_clean$score = rowSums(df_clean[,1:10])

#dichotomization
dich = df_clean
dich[,1:10] = data.frame(lapply(df_clean[,1:10], function (x) as.numeric(x > 2)))
dich$score = rowSums(dich[,1:10])

}

#####
#part 2: CTT-style item analysis
#####
{
  
  #biserial correlations
  biserial_cor = biserial(dich[,SCS_vars],dich[,SCS_vars])
  ggcorrplot(biserial_cor, type = "lower", lab = TRUE)
  ggsave("biserial_cor_mat.pdf",width = 6, height = 6)
  #dichotomous item statistics (percent and N correct, discriminativity)
  dich.distro = rbind(as.character(round(100*unlist(lapply(dich[,SCS_vars], mean)),1)),
                      as.character(as.integer(unlist(lapply(dich[,SCS_vars], sum)))))
  rownames(dich.distro) = c("percent in category 1", "number of cases in category 1")
  
  discrimination = c()
  for (item in 1:10){
    itemname = SCS_vars[item]
    discrimination[itemname] = as.character(round(biserial(rowSums(dich[,-item]),dich[,item]),2))
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
  data_for_eRm = dich[,1:10]
  rasch_model_eRm = RM(data_for_eRm)
  
  
  #approach 2: lme4
  #prepare data for lme4 estimation
  #(item and subject data in long format)
  data_for_lme4 = dich[,1:10]
  data_for_lme4$id = 1:nrow(data_for_lme4)
  data_for_lme4 = melt(data_for_lme4, id.vars = "id")
  rasch_model_lme4 = glmer(value~0+variable+(1|id), data = data_for_lme4,
                        family = binomial)
  smr_lme4 = summary(rasch_model_lme4)
  
  #TODO check: why are estimates different?
  #aproach 3: lavaan
  #modified copy from https://jonathantemplin.com/wp-content/uploads/2022/02/EPSY906_Example05_Binary_IFA-IRT_Models.nb.html
  lavaansyntax = "

    # loadings/discrimination parameters:
    SCS =~ l*Q1 + l*Q2 + l*Q3 + l*Q4 + l*Q5 + l*Q6 + l*Q7 + l*Q8 + l*Q9 + l*Q10
    
    # threshholds use the | operator and start at value 1 after t:
    Q1 | t1; Q2 | t1; Q3 | t1; Q4 | t1; Q5 | t1; Q6 | t1; Q7 | t1; Q8 | t1; Q9 | t1;Q10 | t1;
    
    # factor mean:
    SCS ~ 0;
    
    # factor variance:
    SCS ~~ 1*SCS
    
    "
  data_for_lavaan = dich[,SCS_vars]
  rasch_model_lavaan = sem(model = lavaansyntax, data = data_for_lavaan, ordered = SCS_vars,
                         mimic = "Mplus", estimator = "WLSMV", std.lv = TRUE, parameterization = "theta")
  smr_lavaan = summary(rasch_model_lavaan, fit.measures = TRUE, rsquare = TRUE, standardized = TRUE)
  
  
  #make ICC plot function
  ICC_plot = function(betas){
    df = data.frame(x=seq(-6,6,.01))
    for (i in 1:length(betas)){
      df[[SCS_vars[i]]] = logistic(df$x, betas[i])
    }
    
    df = melt(df, id.vars = "x")
    colnames(df)[2] = "item"
    plt=ggplot(df, aes(x = x, y = value, color = item, label = item)) + geom_line() + 
      theme_clean() + xlab("Person parameter") + ylab("P(item solved)")
    return(directlabels::direct.label(plt, "last.qp"))
  }
  
  #make ICC plots
  betas_eRm = rasch_model_eRm$betapar
  iccplot_eRm=ICC_plot(betas_eRm)
  iccplot_eRm+ggtitle("eRm")

  betas_lme4 = smr_lme4$coefficients[,"Estimate"]
  iccplot_lme4 = ICC_plot(betas_lme4)
  iccplot_lme4+ggtitle("lme4")
  
  betas_lavaan = smr_lavaan$PE[(smr_lavaan$PE$op == "|"),"est"]
  iccplot_lavaan=ICC_plot(betas_lavaan)
  iccplot_lavaan+ggtitle("lavaan")
  
  }

#alternative model: 2PL
{}


#DIF
{}


#reliability, unidimensionality
{}

#measurement invariance
#polytomous IRT model