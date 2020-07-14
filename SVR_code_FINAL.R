# SVR_code_final
# Ola Ozernov-Palchik 
#
# Updated-July 2020
#
# MIT Adult Brain and Cognitive Underpinnings of Dyslexia Study
# SVR paper
#

#### Load and organize ####
Packages <- c("dplyr", "stats", "psych", "ggplot2", "lme4", "lmerTest","Jmisc",
              "gridExtra","olsrr",'relaimpo','BayesFactor','MASS','psych')
lapply(Packages, library, character.only = TRUE)
setwd("~/Dropbox (MIT)/ABCD_analysis") #set working dire

#Load files
d=read.csv("ABCDStudy_2019-09-19.csv")
ev=read.csv("Ev_L_Lang_ROIs.csv")
vwfa=read.csv("vwfa.csv") #vwfa roi for stories > arrow
beh=read.csv("beh_group.csv")

#Clean data
names(beh)[names(beh)=="ids"] <- "PartID"
names(d)[names(d)=="ABCD.ID"] <- "PartID"
names(d)[names(d)=="WRMT.3.Listening.Comprehension..Raw.Score"] <- "LC"
names(d)[names(d)=="GORT.5.Comprehension..Raw.Score"] <- "gortComp"
names(d)[names(d)=="WRMT.3.Word.Attack..Raw.Score"] <- "WA"
names(d)[names(d)=="RAN.RAS.Letters..Raw.Score"] <- "RANL"
names(d)[names(d)=="KBIT.2.Matrices..Standard.Score"] <- "IQ_ss"
names(d)[names(d)=="WRMT.3.Listening.Comprehension..Standard.Score"] <- "LCss"
names(d)[names(d)=="GORT.5.Comprehension..Standard.Score"] <- "gortCompss"
names(d)[names(d)=="WRMT.3.Word.Attack..Standard.Score"] <- "WAss"
names(d)[names(d)=="RAN.RAS.Letters..Raw.Score"] <- "RANL"
names(d)[names(d)=="RAN.RAS.Letters..Standard.Score"] <- "RANLss"
names(d)[names(d)=="KBIT.2.Matrices..Raw.Score"] <- "IQ"
d$Subgroup<-as.factor(d$Subgroup)
d$Sex<-as.factor(d$Sex)



####Behavioral Analysis#####

####run in all subjects

#remove na's
d_lm<-d%>%dplyr::select('PartID','gortComp','Subgroup','Age','Sex','RANL','WA','LC')
d_lm<-na.omit(d_lm)
shapiro.test(d$gortComp)
fit<-lm(gortComp~Age+Sex+RANL+WA+LC, data=d_lm)
ols_vif_tol(fit)
step <- stepAIC(fit, direction = "both",steps = 1000) 
step$anova #
eta_sq(fit)
calc.relimp(fit, type = c("lmg"),
            rela = TRUE)
boot<-boot.relimp(
  fit,
  b = 1000,
  typesel = c("lmg"),
  rank = TRUE,
  diff = TRUE,
  rela = TRUE
)
#
booteval.relimp(boot) # print result
results <- (booteval.relimp(boot, sort = TRUE)) # plot result
plot(results, level = 0.8, names.abbrev = 10,cex=1.5,
     main = "Relative Importance")

####Run only in participants with neuroimaging ####

d<-merge(beh,d_lm,"PartID")
d<-merge(ev,d,"PartID") 
d=merge(vwfa,d,"PartID")

d<-d%>%filter(Subgroup=='Dyslexic'|Subgroup=='Typical') #only DD group

d=d%>%filter(!(PartID %in% c('ABCD_1767','ABCD_1776',
                             'ABCD_1783','ABCD_1763','1702',
                             'ABCD_1798','ABCD_1716','ABCD_1774','ABCD_1753',
                             'ABCD_1732','ABCD_1781'))) #exclude motion

##Seperate by group
d_dys<-d%>%filter(Subgroup=='Dyslexic') #N=18/23
d_typ<-d%>%filter(Subgroup=='Typical') #N=19/21

#Typ
fit_t<-lm(gortComp~Age+Sex+RANL+WA+LC, data=d_typ)
anova(fit_t)
step_t <- stepAIC(fit_t, direction = "both",steps = 1000) 
step_t$anova #:Final Model:gortComp ~ LCss
eta_sq(fit_t)
calc.relimp(fit_t, type = c("lmg"),
            rela = TRUE)
boot_t<-boot.relimp(
  fit_t,
  b = 1000,
  typesel = c("lmg"),
  rank = TRUE,
  diff = TRUE,
  rela = TRUE
)
#
booteval.relimp(boot_t) # print result
results <- (booteval.relimp(boot_t, sort = TRUE)) # plot result
plot(results, level = 0.8, names.abbrev = 10,cex=1.5,
     main = "Relative Importance")

#Dys
fit_d<-lm(gortComp~Age+Sex+RANL+WA+LC, data=d_dys)
step_d <- stepAIC(fit_d, direction = "both",steps = 1000) 
step_d$anova #:Final Model:gortComp ~ LCss
eta_sq(fit_d)
calc.relimp(fit_d, type = c("lmg"),
            rela = TRUE)
boot_d<-boot.relimp(
  fit_d,
  b = 1000,
  typesel = c("lmg"),
  rank = TRUE,
  diff = TRUE,
  rela = TRUE
)

booteval.relimp(boot_d) # print result
results <- (booteval.relimp(boot_d, sort = TRUE)) # plot result
plot(results, level = 0.8, names.abbrev = 10,cex=1.5,
     main = "Relative Importance")

#####Brain Analysis#####

ROI_group <- function(df, dfArr,Range){
  for (i in Range){ 
    fit <- lm(df[,c(i)]~Age+Sex,data=df,na.action=na.exclude)
    df$tmp <- resid(fit)
    norm=shapiro.test(df$tmp)
    if (norm$p.value <=0.05) {
      x <- wilcox.test(df$tmp~dfArr)
    }  else {
      x <- t.test(df$tmp~ dfArr) 
    }
    
    if (x[3] < 0.1) {  
      print(names(df[i]))
      print(x$p.value)
      print(x$statistic)
    }	
  }
}
ROI_group(d, d$Subgroup,2:46) #n.s.

ROI_group2 <- function(df, dfArr,Range){
  for (i in Range){ 
    fit <- lm(df[,c(i)]~Age+Sex,data=df,na.action=na.exclude)
    df$tmp <- resid(fit)
    norm=shapiro.test(df$tmp)
    if (norm$p.value <=0.05) {
      x <- cor.test(df$tmp, dfArr, use="pairwise", 
                    adjust="none",method="spearman", alpha=.05)
    }  else {
      x <- cor.test(df$tmp, dfArr, use="pairwise", 
                    adjust="none",method="pearson", alpha=.05) 
    }
    
    if (x[3] < 0.06) {  
      print(names(df[i]))
      print(x$p.value)
      print(x$estimate)
    }	
  }
}

d1 = d
for (i in 2:37){ 
  fit <- lm(d1[,c(i)]~Age+Sex,data=d1,na.action=na.exclude)
  res <- resid(fit)
  n = paste('res_',names(d1[i]),sep='')
  d1 <- cbind(d1,as.data.frame(res))
  names(d1)[length(names(d1))] <- n
}  


ROI_group2(d, d$LC,2:46) #MD13 -0.29, p=0.08 (L precentral/IFG)
ROI_group2(d, d$RANL,2:46) # 
ROI_group2(d, d$WA,2:46)# *RH_Lang2 -.28 (R MTG/STG), p=0.06, Lang_RH_n1 p=0.06
ROI_group2(d, d$gortComp,2:46) #MD15* -.33, p=0.04 (R SMG/SPL/ANG)

#RANL model
m1<-lm(RANL~Age+Sex+Lang4*Subgroup+Lang1*Subgroup+MD5*Subgroup+vwfa*Subgroup,data=d)
ols_vif_tol(m1)
summary(m1)
step <- stepAIC(m1, direction = "both",na.omit=TRUE)
step$anova #RANL ~ Age + Sex + Subgroup + MD5 + vwfa + Subgroup:MD5

#LC model
m2<-lmp(LC ~ Age + Sex + Subgroup+MD6 +MD1+vwfa+RH_Lang8+RH_Lang1+Lang4+RH_Lang4+
          Lang3+MD6:Subgroup +MD1:Subgroup+vwfa:Subgroup+RH_Lang8:Subgroup+RH_Lang1:Subgroup+
          Lang4:Subgroup+RH_Lang4:Subgroup+Lang3:Subgroup,data=d)
ols_vif_tol(m1)
summary(m2)
step <- stepAIC(m2, direction = "both",na.omit=TRUE)
step$anova
#Final Model:#Final Model:
#LC ~ Age + Sex + Subgroup + MD1 + RH_Lang1 + RH_Lang4 + Lang4 + 
#Subgroup:MD1 + Subgroup:RH_Lang1

m3<-lm(WA~Age+Sex+MD17*Subgroup+ MD5*Subgroup+MD6*Subgroup+Lang5*Subgroup+Lang7*Subgroup+
         RH_Lang5*Subgroup+RH_Lang7*Subgroup,data=d)
ols_vif_tol(m3)
summary(m3)
step <- stepAIC(m3, direction = "both",na.omit=TRUE)
step$anova
#Final Model:
#WA ~ MD17 + Subgroup + MD6 + Lang5 + RH_Lang7 + MD17:Subgroup + 
 # Subgroup:MD6 + Subgroup:Lang5 + Subgroup:RH_Lang7

sqrt(diag(vcov(m1)))
confint(m1)

sqrt(diag(vcov(m2)))
confint(m2)

sqrt(diag(vcov(m3)))
confint(m3)

##### Seperate by group ####

d1$Subgroup<-as.factor(d1$Subgroup)
d_dys<-d%>%filter(Subgroup=='Dyslexic')
d_typ<-d%>%filter(Subgroup=='Typical')

##Typ
ROI_group2(d_typ, d_typ$LC,2:46) # RH_Lang4* OP/ANG -0.59, p=0.008
ROI_group2(d_typ, d_typ$RANL,2:46) #
ROI_group2(d_typ, d_typ$WA,2:46) #
ROI_group2(d_typ, d_typ$gortComp,2:46) #Lang8, RH_Lang4*, MD3*, MD4*, MD16-17*

##models from ROI
m_t1<-lm(RANL~Age+Sex+MD5,data=d_typ)
summary(m_t1)
step <- stepAIC(m_t1)
step$anova #not sig

m_t2<-lmp(LC ~ Age + Sex+MD6 +MD1+vwfa+RH_Lang8+RH_Lang1+Lang4+RH_Lang4+Lang3,data=d_typ)
m_t2<-lmp(LC ~ Age + Sex+MD6 +MD1+vwfa+RH_Lang8+RH_Lang1+Lang4+
            RH_Lang4+Lang3+Lang5+RH_Lang5+Lang7+RH_Lang7,data=d_typ)
summary(m_t2)
step <- stepAIC(m_t2)
step$anova #LC ~ Age + Sex + MD6 + MD1 + vwfa + RH_Lang8 + Lang4 + RH_Lang4 + Lang3

m_t3<-lm(WA~Age+Sex+MD6+Lang5+RH_Lang7,data=d_typ) 
summary(m_t3)
step <- stepAIC(m_t3)
step$anova #WA ~ Sex + MD6 + Lang5 + RH_Lang7

#predicting reading comp
fit_t<-lm(gortComp~Age+RH_Lang4+Lang8+MD3+MD4+MD16,data=d_typ) #RANL ~ Age + Sex + Subgroup + MD5 + vwfa + Subgroup:MD5
step <- stepAIC(fit_t, direction = "both",steps = 1000)
step$anova
cor.test(d_typ$vwfa,d_typ$RANL,use="pairwise", method="pearson", alpha=.05)

##Dys
ROI_group2(d_dys, d_dys$LC,2:46) #all positive: RHLang1*, RHLang3, RHLang5*,MD9*,Lang_RH_n1*,Lang_RH_n5*
ROI_group2(d_dys, d_dys$RANL,2:46) #all neg: Lang6, RHLang4*,RHLang7*,Lang_RH_n4, MD4*,MD8*,MD11,MD17
ROI_group2(d_dys, d_dys$WA,2:46) #all neg Lang1*, Lang2*, Lang5*, RHLang1*, RHLang5*, MD1*,MD9*,MD10*,Lang_n1,Lang_RH_n1
ROI_group2(d_dys, d_dys$gortComp,2:46) #pos: RHLang3, MD17

anova(lm(RANL~Age+Sex+vwfa+Lang5+RH_Lang7+MD4+MD5,data=d_dys))
fit_dys<-lm(gortComp~Age+Sex+WA+LC+RANL+RH_Lang3+MD17,data=d_dys)

step <- stepAIC(fit_dys, direction = "both",steps = 1000)
step$anova
summary(lm(gortComp ~ Age+Sex+RH_Lang3,data=d_dys))
##ROI
m_d1<-lm(RANL~Age+Sex+MD5,data=d_dys) #all sig!
summary(m_d1)
step <- stepAIC(m_d1)
step$anova

m_d2<-lmp(LC ~ Age + Sex+MD6 +MD1+vwfa+RH_Lang8+RH_Lang1+Lang4+RH_Lang4+
            Lang3,data=d_dys)
m_d2<-lmp(LC ~ Age + Sex+MD6 +MD1+vwfa+RH_Lang8+RH_Lang1+Lang4+
            RH_Lang4+Lang3+Lang5+RH_Lang5+Lang7+RH_Lang7,data=d_dys)

summary(m_d2)
step <- stepAIC(m_d2)
step$anova #LC ~ RH_Lang1 + Lang4
sqrt(diag(vcov(m1)))


m_d3<-lm(WA~Age+Sex+MD6+Lang5+RH_Lang7,data=d_dys) 
summary(m_d3)
step <- stepAIC(m_d3)
step$anova #WA ~ MD6 + Lang5 + RH_Lang7



########ggplot

ggplot(d1, aes(x = LC, y = res_MD1)) +
  stat_smooth (
    method = "glm",
    formula = y ~ x,
    # colour = "black",
    size = 1
  ) +
  geom_point(aes(shape = as.factor(Subgroup)), size = 3) +
  scale_shape_manual(values = c(1, 17)) +
  labs(x = "LC", y = "MD17") +
  theme(
    axis.title = element_text(family = "Trebuchet MS", size = 20),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(
      size = 10,
      face = 'bold' ),
    legend.title = element_blank(),
    #legend.position = c(0.85, 0.15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )

ggplot(d1, aes(x=LC,y=res_MD1,shape = Subgroup, linetype=Subgroup)) +
  geom_point(size = 5) +
  geom_smooth (method="lm",col="red")+
  scale_shape_manual(values = c(1, 17)) +
  labs(x = "LC", y = "res_MD1") +
  theme(
    axis.title = element_text(family = "Trebuchet MS", size = 20),
    legend.position ="none",
    legend.text = element_blank(),
    legend.title = element_blank(),
    panel.background = element_rect(fill = "transparent", colour = NA), 
    plot.background = element_rect(fill = "transparent", colour = NA),
    axis.line = element_line(colour = "black"),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15))

d1_dys<-d1%>%filter(Subgroup=='Dyslexic')
d1_typ<-d1%>%filter(Subgroup=='Typical')

ggplot(d1, aes(x = LC, y = res_MD17)) +
  stat_smooth (
    method = "glm",
    formula = y ~ x,
    # colour = "black",
    size = 1
  ) +
  geom_point(size = 5) +
  scale_shape_manual(values = c(1, 17)) +
  labs(x = "LC", y = "MD17") +
  theme(
    axis.title = element_text(family = "Trebuchet MS", size = 20),
    legend.key.size = unit(1.5, "cm"),
    legend.text = element_text(
      size = 10,
      face = 'bold' ),
    legend.title = element_blank(),
    #legend.position = c(0.85, 0.15),
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15)
  )
