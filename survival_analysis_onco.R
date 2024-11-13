###
# Recurrence of Clear Cell Renal Cell Carcinoma after Treatment
# author: Rafa Carretero
# date: 11 november 2024
###
# 
#
#######################################
#                                     #
# LOADING LIBRARIES                   #
#                                     #
#######################################
#
#Clear existing data and graphics
rm(list=ls())
graphics.off()

library(MASS)
library(glmnet)
library(ggplot2)
library(gtsummary)
library(dplyr)
library(Hmisc)
library(scales)
library(viridis)
library(corrplot)
library(ggsignif)
# librerías para el 
# análisis de supervivencia:
library(survival)
library(ggsurvfit)
library(survminer)

#######################################
#                                     #
# DATA RECOVERY                       #
#                                     #
#######################################
#
setwd("~/Python_Projects/ONCO_Genes/")
pacientes <- read.csv("data_onco.csv", header=TRUE, sep=",")

# Data pre-processing
pacientes$AXL <- as.numeric(gsub(",", ".", gsub("\\.", "", pacientes$AXL)))
pacientes$BAP1 <- as.numeric(gsub(",", ".", gsub("\\.", "", pacientes$BAP1)))
pacientes$CA9 <- as.numeric(gsub(",", ".", gsub("\\.", "", pacientes$CA9)))

pacientes$Tipo_Histologico[pacientes$Tipo_HistologicoRCC_0ccRCC_1Otro == 1] <- "Otro"
pacientes$Tipo_Histologico[pacientes$Tipo_HistologicoRCC_0ccRCC_1Otro == 0] <- "Ca. células claras"

pacientes$Sexo_H0_M1[pacientes$Sexo_H0_M1 == 1] <- "Mujer"
pacientes$Sexo_H0_M1[pacientes$Sexo_H0_M1 == 0] <- "Hombre"

pacientes$Mejor_Respuesta[pacientes$MejorRespuesta_1erTKI == "RC"] <- "Sí"
pacientes$Mejor_Respuesta[pacientes$MejorRespuesta_1erTKI == "RP"] <- "Sí"
pacientes$Mejor_Respuesta[pacientes$MejorRespuesta_1erTKI == "EE"] <- "Sí"
pacientes$Mejor_Respuesta[pacientes$MejorRespuesta_1erTKI == "PE"] <- "No"

pacientes$Progresion_Enfermedad[pacientes$PE_1erTKI_0No_1Si == 1] <- "Sí"
pacientes$Progresion_Enfermedad[pacientes$PE_1erTKI_0No_1Si == 0] <- "No"
pacientes$Progresion_Enfermedad <- as.factor(pacientes$Progresion_Enfermedad)

#######################################
#                                     #
# TABLA DESCRIPTIVA (I)               #
#                                     #
#######################################
#
pacientes %>%
  tbl_summary()%>%
  bold_labels()

pacientes %>%
  select(sexo,
         #Sexo_H0_M1,
         Age_1erTKI,
         GrupoPronostico_MSKCC_en1lineaTKI, 
         #GrupoPronostico_MSKCC_0Bueno_1Intermedio_2Malo,
         X1erTKI_cual,
         #MejorRespuesta_1erTKI,
         Mejor_Respuesta,
         Tipo_Histologico,
         nefrectomia,
         Progresion_Enfermedad,
         PFS_1erTKI_days,
         OS_days,
         EXITUS_si_no) %>%
  tbl_summary(
    #by=Progresion_Enfermedad,
    label = list(sexo ~ "Sexo",
                 Age_1erTKI ~ "Edad",
                 GrupoPronostico_MSKCC_en1lineaTKI ~ "Grupo pronóstico",
                 X1erTKI_cual ~ "Tratamiento Primera Línea", 
                 Tipo_Histologico ~ "Tipo Histológico", 
                 Mejor_Respuesta ~ "Respuesta tras Primer Tratamiento", 
                 nefrectomia ~ "Nefrectomía",
                 PFS_1erTKI_days ~ "Tiempo Libre de Enfermedad",
                 Progresion_Enfermedad ~ "Progresión de la Enfermedad",
                 OS_days ~ "Tiempo de Supervivencia Global",
                 EXITUS_si_no ~ "Fallecimiento"),
  ) %>%
  #add_overall() %>%
  modify_footnote(
    all_stat_cols() ~ "Medidas expresadas en Mediana (RIQ) o Frecuencia (%)"
  ) %>%
  modify_caption("**Tabla 1. Características generales de la cohorte de estudio**") %>%
  modify_header(label ~ "**Variable**") %>%
  #modify_spanning_header(c("stat_1", "stat_2") ~ "**Progresión de la Enfermedad**") %>%
  bold_labels()
#
#######################################
#                                     #
# TABLA DESCRIPTIVA (II)              #
#                                     #
#######################################
#
pacientes %>%
  select(sexo,
         #Sexo_H0_M1,
         Age_1erTKI,
         GrupoPronostico_MSKCC_en1lineaTKI, 
         #GrupoPronostico_MSKCC_0Bueno_1Intermedio_2Malo,
         X1erTKI_cual,
         #MejorRespuesta_1erTKI,
         Mejor_Respuesta,
         Tipo_Histologico,
         nefrectomia,
         Progresion_Enfermedad,
         PFS_1erTKI_days,
         OS_days,
         EXITUS_si_no) %>%
  tbl_summary(
    by=Progresion_Enfermedad,
    label = list(sexo ~ "Sexo",
                 Age_1erTKI ~ "Edad",
                 GrupoPronostico_MSKCC_en1lineaTKI ~ "Grupo pronóstico",
                 X1erTKI_cual ~ "Tratamiento Primera Línea", 
                 Tipo_Histologico ~ "Tipo Histológico", 
                 Mejor_Respuesta ~ "Respuesta tras Primer Tratamiento", 
                 nefrectomia ~ "Nefrectomía",
                 PFS_1erTKI_days ~ "Tiempo Libre de Enfermedad",
                 #Progresion_Enfermedad ~ "Progresión de la Enfermedad",,
                 OS_days ~ "Tiempo de Supervivencia Global",
                 EXITUS_si_no ~ "Fallecimiento"),
  ) %>%
  add_overall() %>%
  modify_footnote(
    all_stat_cols() ~ "Medidas expresadas en Mediana (RIQ) o Frecuencia (%)"
  ) %>%
  modify_caption("**Tabla 2. Características generales de la cohorte de estudio**") %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Progresión de la Enfermedad**") %>%
  bold_labels()
#
#######################################
#                                     #
# EXPRESION ARNm - AXL                #
#                                     #
#######################################
#
pacientes$AXL
boxplot(pacientes$AXL, pacientes$Progresion_Enfermedad)

#Figura 1
#########
ggplot(pacientes, aes(x=AXL)) + 
  geom_histogram(aes(y=after_stat(density)),      # Histogram with density instead of count on y-axis
                 binwidth=.2,
                 colour=rgb(0.2,0.8,0.5,0.5), fill=rgb(0.2,0.8,0.5,0.5)) +
  geom_density(alpha=.2, fill="#2D8B3F") + # Overlay with transparent density plot
  theme_bw() +
  ggtitle("Expresión de ARNm para AXL", subtitle = "Estudio descriptivo (mediana)") +
  geom_vline(aes(xintercept=median(AXL, na.rm=T)),   # Ignore NA values for median
             color="red", linetype="dashed", linewidth=0.5) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.title=element_blank())
ggsave("Figuras/Figura_1.png", height = 6, width = 8)

#Figura 2
#########
ggplot(pacientes, aes(x = 0, y = AXL)) + 
  geom_boxplot(colour="#09622A", fill=rgb(0.2,0.8,0.5,0.5) )+
  geom_jitter() +
  theme_bw() +
  ggtitle("Expresión de ARNm para AXL", subtitle = "Estudio descriptivo (boxplot)") +
  #scale_fill_brewer(palette="Paired") +
  labs(x='', y='Expresión ARNm de AXL') +
  theme(legend.title=element_blank())
ggsave("Figuras/Figura_2.png", height = 6, width = 8)

#
#######################################
#                                     #
# EXPRESION ARNm - BAP1               #
#                                     #
#######################################
#
pacientes$BAP1
boxplot(pacientes$BAP1, pacientes$Progresion_Enfermedad)

#Figura 3
#########
ggplot(pacientes, aes(x=BAP1)) + 
  geom_histogram(aes(y=after_stat(density)),      # Histogram with density instead of count on y-axis
                 binwidth=.3,
                 colour="#8DBBDC", fill="#99C5E3") +
  geom_density(alpha=.2, fill="#B9DDF1") + # Overlay with transparent density plot
  theme_bw() +
  ggtitle("Expresión de ARNm para BAP1", subtitle = "Estudio descriptivo (mediana)") +
  geom_vline(aes(xintercept=median(AXL, na.rm=T)),   # Ignore NA values for median
             color="red", linetype="dashed", linewidth=0.5) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.title=element_blank())
ggsave("Figuras/Figura_3.png", height = 6, width = 8)

#Figura 4
#########
ggplot(pacientes, aes(x = 0, y = BAP1)) + 
  geom_boxplot(colour="black", fill="#8DBBDC")+
  geom_jitter() +
  theme_bw() +
  ggtitle("Expresión de ARNm para BAP1", subtitle = "Estudio descriptivo (boxplot)") +
  #scale_fill_brewer(palette="Paired") +
  labs(x='', y='Expresión ARNm de BAP1') +
  scale_y_continuous(limits = c(0, 3.2)) +
  theme(legend.title=element_blank())
ggsave("Figuras/Figura_4.png", height = 6, width = 8)

#
#######################################
#                                     #
# EXPRESION ARNm - CA9                #
#                                     #
#######################################
#
pacientes$CA9
boxplot(pacientes$CA9, pacientes$Progresion_Enfermedad)

#Figura 5
#########
ggplot(pacientes, aes(x=CA9)) + 
  geom_histogram(aes(y=after_stat(density)),      # Histogram with density instead of count on y-axis
                 binwidth=.3,
                 colour="#BB173A", fill="#F26553") +
  geom_density(alpha=.2, fill="#F46C58") + # Overlay with transparent density plot
  theme_bw() +
  ggtitle("Expresión de ARNm para CA9", subtitle = "Estudio descriptivo (mediana)") +
  geom_vline(aes(xintercept=median(AXL, na.rm=T)),   # Ignore NA values for median
             color="red", linetype="dashed", linewidth=0.5) +
  scale_fill_brewer(palette="Paired") +
  theme(legend.title=element_blank())
ggsave("Figuras/Figura_5.png", height = 6, width = 8)

#Figura 6
#########
ggplot(pacientes, aes(x = 0, y = CA9)) + 
  geom_boxplot(colour="black", fill="#F26553")+
  geom_jitter() +
  theme_bw() +
  ggtitle("Expresión de ARNm para CA9", subtitle = "Estudio descriptivo (boxplot)") +
  #scale_fill_brewer(palette="Paired") +
  labs(x='', y='Expresión ARNm de CA9') +
  scale_y_continuous(limits = c(0, 4)) +
  theme(legend.title=element_blank())
ggsave("Figuras/Figura_6.png", height = 6, width = 8)

#
#######################################
#                                     #
# CATEGORIAS EXPRESION ARNm           #
#                                     #
#######################################
#
pacientes$AXL_cat <- cut(pacientes$AXL,
                     breaks=c(-1, 1, 3),
                     labels=c('<Mediana', '>Mediana'))
pacientes$BAP1_cat <- cut(pacientes$BAP1,
                          breaks=c(-1, 1, 5.4),
                          labels=c('<Mediana', '>Mediana'))
pacientes$CA9_cat <- cut(pacientes$CA9,
                         breaks=c(-1, 1, 5.5),
                         labels=c('<Mediana', '>Mediana'))

#
#######################################
#                                     #
# TABLA DESCRIPTIVA (III)             #
#                                     #
#######################################
#
pacientes %>%
  select(sexo,
         #Sexo_H0_M1,
         Age_1erTKI,
         GrupoPronostico_MSKCC_en1lineaTKI, 
         #GrupoPronostico_MSKCC_0Bueno_1Intermedio_2Malo,
         X1erTKI_cual,
         #MejorRespuesta_1erTKI,
         Mejor_Respuesta,
         Tipo_Histologico,
         AXL_cat,
         BAP1_cat,
         CA9_cat,
         nefrectomia,
         Progresion_Enfermedad,
         PFS_1erTKI_days,
         OS_days,
         EXITUS_si_no) %>%
  tbl_summary(
    by=Progresion_Enfermedad,
    label = list(sexo ~ "Sexo",
                 Age_1erTKI ~ "Edad",
                 GrupoPronostico_MSKCC_en1lineaTKI ~ "Grupo pronóstico",
                 X1erTKI_cual ~ "Tratamiento Primera Línea", 
                 Tipo_Histologico ~ "Tipo Histológico", 
                 AXL_cat ~ "Expresión AXL",
                 BAP1_cat ~ "Expresión BAP1",
                 CA9_cat ~ "Expresión CA9",
                 Mejor_Respuesta ~ "Respuesta tras Primer Tratamiento", 
                 nefrectomia ~ "Nefrectomía",
                 PFS_1erTKI_days ~ "Tiempo Libre de Enfermedad",
                 #Progresion_Enfermedad ~ "Progresión de la Enfermedad",,
                 OS_days ~ "Tiempo de Supervivencia Global",
                 EXITUS_si_no ~ "Fallecimiento"),
  ) %>%
  add_p(pvalue_fun = ~ style_pvalue(.x, digits = 2)) %>%
  add_overall() %>%
  modify_footnote(
    all_stat_cols() ~ "Medidas expresadas en Mediana (RIQ) o Frecuencia (%)"
  ) %>%
  modify_caption("**Tabla 3. Características generales de la cohorte de estudio (con expresión ARNm)**") %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(c("stat_1", "stat_2") ~ "**Progresión de la Enfermedad**") %>%
  bold_labels() %>%
  italicize_levels() %>%
  bold_p(t = 0.05)
#
#######################################
#                                     #
# REGRESION LOGISTICA                 #
#                                     #
#######################################
# En grupo pronóstico agrupo intermedio/malo, 
# porque hay sólo dos pacientes en el grupo "malo"
pacientes$GrupoPronostico_MSKCC_en1lineaTKI[pacientes$GrupoPronostico_MSKCC_en1lineaTKI == "Intermedio"] <- "Intermedio/Malo"
pacientes$GrupoPronostico_MSKCC_en1lineaTKI[pacientes$GrupoPronostico_MSKCC_en1lineaTKI == "Malo"] <- "Intermedio/Malo"

pacientes$sexo = relevel(as.factor(pacientes$sexo), ref="Mujer")
#pacientes$AXL_cat = relevel(as.factor(pacientes$AXL_cat), ref=">Mediana")

# usamos la variable: pacientes$PE_1erTKI_0No_1Si
model1 <- glm(PE_1erTKI_0No_1Si ~ sexo + 
                Age_1erTKI + GrupoPronostico_MSKCC_en1lineaTKI + 
              X1erTKI_cual + Mejor_Respuesta + Tipo_Histologico + 
              AXL_cat +
              BAP1_cat + CA9_cat +
              #AXL + BAP1 + CA9 +
              nefrectomia, pacientes, family = binomial)

tbl_regression(model1, 
               exponentiate = TRUE,
               label = list(sexo	 ~ "Sexo",
                            #AXL_cat ~ "AXL_cat",
                            Age_1erTKI ~ "Edad",
                            GrupoPronostico_MSKCC_en1lineaTKI ~ "Grupo pronóstico",
                            X1erTKI_cual ~ "Tratamiento Primera Línea", 
                            Tipo_Histologico ~ "Tipo Histológico", 
                            Mejor_Respuesta ~ "Respuesta tras Primer Tratamiento", 
                            nefrectomia ~ "Nefrectomía"
                            ),
               ) %>%
  #add_nevent() %>% # add number of events of the outcome
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Tabla 4. Análisis multivariado con regresión logística**") %>%
  modify_header(label ~ "**Variable**") %>%
  italicize_levels()

#
#######################################
#                                     #
# REGRESION LOGISTICA - STEPWISE      #
#                                     #
#######################################
#
step.model <- stepAIC(model1, direction = "both", 
                      trace = FALSE)
summary(step.model)
#glm(formula = PE_1erTKI_0No_1Si ~ sexo + GrupoPronostico_MSKCC_en1lineaTKI + 
#      X1erTKI_cual + AXL_cat + BAP1_cat, family = binomial, data = pacientes)


tbl_regression(step.model, 
               exponentiate = TRUE,
               label = list(sexo	 ~ "Sexo",
                            GrupoPronostico_MSKCC_en1lineaTKI ~ "Grupo pronóstico",
                            X1erTKI_cual ~ "Tratamiento Primera Línea"
               ),
) %>%
  #add_nevent() %>% # add number of events of the outcome
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Tabla 5. Análisis multivariado con regresión logística (stepwise)**") %>%
  modify_header(label ~ "**Variable**") %>%
  italicize_levels()

#
#######################################
#                                     #
# REGRESION LOGISTICA - LASSO         #
#                                     #
#######################################
#
# Prepare predictor matrix and outcome variable
x <- model.matrix(PE_1erTKI_0No_1Si ~ sexo + 
                    Age_1erTKI + GrupoPronostico_MSKCC_en1lineaTKI + 
                    X1erTKI_cual + Mejor_Respuesta + Tipo_Histologico + 
                    AXL_cat +
                    BAP1_cat + CA9_cat +nefrectomia, data = pacientes)[, -1]
y <- as.numeric(pacientes$PE_1erTKI_0No_1Si)  # Ensure binary outcome is numeric (0,1)

# Perform Lasso logistic regression
fit.lasso <- glmnet(x, y, family="binomial", alpha=1)
lasso_model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

# Best lambda value
best_lambda_1 <- lasso_model$lambda.min
best_lambda_2 <- lasso_model$lambda.1se

best_lambda_1
best_lambda_2

# Coefficients of the Lasso model with best lambda
coef(lasso_model, s = best_lambda_1)
coef(lasso_model, s = best_lambda_2)

plot(lasso_model)
plot(fit.lasso)

#
#######################################
#                                     #
# CURVA SUPERVIVENCIA (I)             #
#                                     #
#######################################
#
# la variable "recaida" (progresión de la enfermedad)
# tiene que volver a ser
# número, porque si no, el análisis de supervivencia no sale
# usamos la variable: pacientes$PE_1erTKI_0No_1Si

km_fit_1 <- survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ 1, data=pacientes)

curva_1 <- ggsurvplot(km_fit_1, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = TRUE,
                      #pval=TRUE,
                      surv.median.line = "hv", # Specify median survival
                      ggtheme = theme_bw(),
                      legend = c("bottom"),
                      title = "Curva de Kaplan-Meier - Tiempo libre de enfermedad",
                      xlab = "Tiempo (días)",
                      ylab = "Libre de Enfermedad (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      censor = FALSE
)
curva_1$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  scale_colour_manual(values = c("steelblue","red","red","red")) +
  labs(title = "Curva de Kaplan-Meier",
       subtitle="Tiempo libre de enfermedad")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))

ggsave("Figuras/Figura_7.png", height = 6, width = 8)

tbl_survfit(
  survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ 1, data=pacientes),
  probs = 0.5,
  label_header = "**Mediana de Surpervivencia**"
) %>%
  #bold_labels() %>%
  modify_header(label ~ "**Supervivencia**") %>%
  italicize_levels()

#
#######################################
#                                     #
# CURVA SUPERVIVENCIA (II)            #
#                                     #
#######################################
#
# la variable "recaida" (progresión de la enfermedad)
# tiene que volver a ser
# número, porque si no, el análisis de supervivencia no sale
# usamos la variable: pacientes$PE_1erTKI_0No_1Si

km_fit_1 <- survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ AXL_cat, data=pacientes)

curva_1 <- ggsurvplot(km_fit_1, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = TRUE,
                      #pval=TRUE,
                      pval="log-rank test \n p=0.06",
                      surv.median.line = "hv", # Specify median survival
                      palette = c("#E7B800","#2E9FDF"),# custom color palettes
                      ggtheme = theme_bw(),
                      legend = c("bottom"),
                      title = "Estimador de Kaplan-Meier - Tiempo libre de enfermedad",
                      xlab = "Tiempo (días)",
                      ylab = "Libre de Enfermedad (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      legend.title="", # elimina "strata"
                      risk.table.col = "strata",# Risk table color by groups
                      #legend.labs = c("Baja expresión AXL", "Alta expresión AXL"),    # Change legend labels
                      censor = FALSE,
                      #break.time.by=730.5 #cada 2 años, contando el bisiesto
)
curva_1$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  #scale_colour_manual(values = c("steelblue","red","red","red")) +
  labs(title = "Estimador de Kaplan-Meier - Expresión de AXL",
       subtitle="Tiempo libre de enfermedad")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))

ggsave("Figuras/Figura_8.png", height = 6, width = 8)

tbl_survfit(
  survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ AXL_cat, data=pacientes),
  probs = 0.5,
  label_header = "**Mediana de Surpervivencia**"
) %>%
  #bold_labels() %>%
  modify_header(label ~ "**Supervivencia**") %>%
  italicize_levels()

#
#######################################
#                                     #
# CURVA SUPERVIVENCIA (III)           #
#                                     #
#######################################
#
# la variable "recaida" (progresión de la enfermedad)
# tiene que volver a ser
# número, porque si no, el análisis de supervivencia no sale
# usamos la variable: pacientes$PE_1erTKI_0No_1Si

km_fit_1 <- survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ BAP1_cat, data=pacientes)

curva_1 <- ggsurvplot(km_fit_1, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = TRUE,
                      #pval=TRUE,
                      pval="log-rank test \n p=0.16",
                      surv.median.line = "hv", # Specify median survival
                      palette = c("#E7B800", "#AF4623"),# custom color palettes
                      ggtheme = theme_bw(),
                      legend = c("bottom"),
                      title = "Estimador de Kaplan-Meier - Tiempo libre de enfermedad",
                      xlab = "Tiempo (días)",
                      ylab = "Libre de Enfermedad (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      legend.title="", # elimina "strata"
                      risk.table.col = "strata",# Risk table color by groups
                      legend.labs = c("Baja expresión BAP1", "Alta expresión BAP1"),    # Change legend labels
                      censor = FALSE
)
curva_1$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  #scale_colour_manual(values = c("steelblue","red","red","red")) +
  labs(title = "Estimador de Kaplan-Meier - Expresión de BAP1",
       subtitle="Tiempo libre de enfermedad")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))

ggsave("Figuras/Figura_9.png", height = 6, width = 8)

#
#######################################
#                                     #
# CURVA SUPERVIVENCIA (IV)            #
#                                     #
#######################################
#
# la variable "recaida" (progresión de la enfermedad)
# tiene que volver a ser
# número, porque si no, el análisis de supervivencia no sale
# usamos la variable: pacientes$PE_1erTKI_0No_1Si

km_fit_1 <- survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ CA9_cat, data=pacientes)

curva_1 <- ggsurvplot(km_fit_1, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = TRUE,
                      #pval=TRUE,
                      pval="log-rank test \n p=0.018",
                      surv.median.line = "hv", # Specify median survival
                      palette = c("#E7B800", "#1A7332"),# custom color palettes
                      ggtheme = theme_bw(),
                      legend = c("bottom"),
                      title = "Estimador de Kaplan-Meier - Tiempo libre de enfermedad",
                      xlab = "Tiempo (días)",
                      ylab = "Libre de Enfermedad (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      legend.title="", # elimina "strata"
                      risk.table.col = "strata",# Risk table color by groups
                      legend.labs = c("Baja expresión CA9", "Alta expresión CA9"),    # Change legend labels
                      censor = FALSE
)
curva_1$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  #scale_colour_manual(values = c("steelblue","red","red","red")) +
  labs(title = "Estimador de Kaplan-Meier - Expresión de CA9",
       subtitle="Tiempo libre de enfermedad")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))

ggsave("Figuras/Figura_10.png", height = 6, width = 8)

#
#######################################
#                                     #
# CURVA SUPERVIVENCIA (V)             #
#                                     #
#######################################
#
# la variable "recaida" (progresión de la enfermedad)
# tiene que volver a ser
# número, porque si no, el análisis de supervivencia no sale
# usamos la variable: pacientes$PE_1erTKI_0No_1Si

km_fit_1 <- survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ BestResp_RCRP0vsPEEE1, data=pacientes)

curva_1 <- ggsurvplot(km_fit_1, 
                      data = pacientes, 
                      risk.table = TRUE, 
                      conf.int = TRUE,
                      #pval=TRUE,
                      pval="log-rank test \n p=0.001",
                      surv.median.line = "hv", # Specify median survival
                      palette = c("#26456E", "#E7B800"),# custom color palettes
                      ggtheme = theme_bw(),
                      legend = c("bottom"),
                      title = "Estimador de Kaplan-Meier - Tiempo libre de enfermedad",
                      xlab = "Tiempo (días)",
                      ylab = "Libre de Enfermedad (probabilidad)",
                      risk.table.title="Pacientes en riesgo",
                      legend.title="", # elimina "strata"
                      risk.table.col = "strata",# Risk table color by groups
                      legend.labs = c("Respuesta Parcial o Completa", "Progresión"),    # Change legend labels
                      censor = FALSE
)
curva_1$plot + 
  theme_bw() +
  scale_linetype_manual(values = c("solid","dashed", "solid", "dotted", "dotdash")) +
  #scale_colour_manual(values = c("steelblue","red","red","red")) +
  labs(title = "Estimador de Kaplan-Meier - Respuesta al Tratamiento",
       subtitle="Tiempo libre de enfermedad")+  
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position="bottom") +
  theme(legend.title = element_blank()) +
  theme(panel.background = element_rect(colour = "black"),
        axis.text=element_text(size=8),  #tamaño de las fechas
        axis.title.x = element_text(vjust=-0.2),
        axis.title.y = element_text(vjust=+0.6),
        axis.title=element_text(size=10,face="bold"), #tamaño de los títulos de los ejes
        plot.title = element_text(size = 11, hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey40"),
        plot.caption = element_text(size = 7.5, color = "grey40"))

ggsave("Figuras/Figura_11.png", height = 6, width = 8)


tbl_survfit(
  survfit(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si) ~ BestResp_RCRP0vsPEEE1, data=pacientes),
  probs = 0.5,
  label_header = "**Mediana de Surpervivencia**"
) %>%
  #bold_labels() %>%
  modify_header(label ~ "**Supervivencia**") %>%
  italicize_levels()


#
#######################################
#                                     #
# REGRESION DE COX                    #
#                                     #
#######################################
#
# Volvemos a poner la referencia de AXL en su sitio:
#pacientes$AXL_cat = relevel(as.factor(pacientes$AXL_cat), ref="<Mediana")

cox.clasico <- coxph(Surv(PFS_1erTKI_days, PE_1erTKI_0No_1Si)~ 
                       #sexo + 
                       #Age_1erTKI + GrupoPronostico_MSKCC_en1lineaTKI + 
                       #Tipo_Histologico + 
                       Mejor_Respuesta +
                       AXL_cat +
                       BAP1_cat + CA9_cat,
                       #+nefrectomia, 
                       data = pacientes,na.action=na.exclude)
summary(cox.clasico)




tabla_cox <- tbl_regression(cox.clasico, exponentiate = TRUE) %>%
  #add_n() %>%
  modify_caption("**Tabla 6. Hazard Ratios Ajustados**") %>%
  modify_header(label = "**Variable**") %>% # update the column header
  bold_labels() %>%
  bold_p(t = 0.05) %>%
  italicize_levels()
tabla_cox


############
############
step.model.cox <- stepAIC(cox.clasico, direction = "both")
summary(step.model.cox)


tbl_regression(step.model.cox, 
               exponentiate = TRUE,
               label = list(sexo	 ~ "Sexo",
                            GrupoPronostico_MSKCC_en1lineaTKI ~ "Grupo pronóstico",
                            X1erTKI_cual ~ "Tratamiento Primera Línea"
               ),
) %>%
  #add_nevent() %>% # add number of events of the outcome
  bold_p(t = 0.05) %>%
  bold_labels() %>%
  modify_caption("**Tabla 6. Hazard Ratios Ajustados**") %>%
  modify_header(label ~ "**Variable**") %>%
  italicize_levels()