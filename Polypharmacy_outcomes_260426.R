# =============================================================================
# Polypharmacy Clinical Outcomes Analysis
# -----------------------------------------------------------------------------
# Reproducibility notes for publication:
# - This script expects original source datasets that are not distributed here.
# - Configure input file locations using environment variables:
#     POLYPHARMACY_MAIN_COHORT_SAV
#     POLYPHARMACY_MEDICATION_SAV
#     POLYPHARMACY_POST_AF_CATEGORY_SAV
# - If not set, the script falls back to legacy local filenames.
# =============================================================================

# Load required packages
library(haven)
library(cohorttools)
library(dplyr)    # ensure dplyr functions are registered last

main_cohort_path <- Sys.getenv(
  "POLYPHARMACY_MAIN_COHORT_SAV",
  unset = "IncidentCohort280924.sav"
)
medication_path <- Sys.getenv(
  "POLYPHARMACY_MEDICATION_SAV",
  unset = "laakkeet120pvEnnenAFkohorttiintuloaKaikki.sav"
)
post_af_category_path <- Sys.getenv(
  "POLYPHARMACY_POST_AF_CATEGORY_SAV",
  unset = "post_af_polypharmacy_categories.sav"
)

#-----------------------------
# Load main cohort dataset
#-----------------------------
cohort_df <- read_sav(main_cohort_path) %>%
  as.data.frame()

# Quick preview
head(cohort_df)

#-----------------------------
# Load medication dataset
#-----------------------------
medication_df <- read_sav(medication_path)

# Quick preview
head(medication_df)

#-----------------------------
# Prepare medication data
#-----------------------------

# Ensure ATC codes are strings and create ATC level-4 variable
medication_df <- medication_df %>%
  dplyr::mutate(
    atc = as.character(atc),
    ATC_4 = substr(atc, 1, 4)
  )

# Keep only relevant columns and remove duplicates (unique ATC per SID)
medication_unique_df <- data %>%
  dplyr::select(SID, atc, ATC_4) %>%
  dplyr::distinct()

# Summarize number of unique drugs per patient
medication_summary_df <- medication_unique_df %>%
  dplyr::group_by(SID) %>%
  dplyr::summarise(
    Unique_Drugs_ATC = dplyr::n_distinct(atc),
    Unique_Drugs_ATC4 = dplyr::n_distinct(ATC_4),
    .groups = "drop"
  )

#-----------------------------
# Merge summaries into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(medication_summary_df, by = "SID") %>%
  dplyr::mutate(
    Unique_Drugs_ATC = ifelse(is.na(Unique_Drugs_ATC), 0, Unique_Drugs_ATC),
    Unique_Drugs_ATC4 = ifelse(is.na(Unique_Drugs_ATC4), 0, Unique_Drugs_ATC4)
  )

#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC, Unique_Drugs_ATC4) %>%
  head()


library(dplyr)

# Categorize Unique_Drugs_ATC into 4 groups
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = dplyr::case_when(
      Unique_Drugs_ATC >= 0  & Unique_Drugs_ATC <= 2  ~ "0–2 drugs",
      Unique_Drugs_ATC >= 3  & Unique_Drugs_ATC <= 4  ~ "3–4 drugs",
      Unique_Drugs_ATC >= 5  & Unique_Drugs_ATC <= 9  ~ "5–9 drugs",
      Unique_Drugs_ATC >= 10                          ~ "10+ drugs",
      TRUE                                            ~ NA_character_   # fallback
    )
  )

# Optional: make it an ordered factor
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = factor(
      Unique_Drugs_ATC_Category,
      levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
      ordered = TRUE
    )
  )

# Quick check
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_Category,
  levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)

with(cohort_df,range(cal.yr(DateISorLoppuOrDeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","IS")))
exit_status_factor<-with(cohort_df,factor(ifelse(ISaftercohortall==1,"IS","EOF"),levels=c("SOF","EOF","IS")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateISorLoppuOrDeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)


boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))



library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Unique_Drugs_ATC_Category, data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTablepolypharmacy_IS.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="IS",lex.dur) ~ Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+
                sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+ALKO+
                dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiom <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiom <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiom <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)


sink(file="IRRTulosUnique_Drugs_ATC_Category_IS.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()

apu<-unclass(lexis_df3$low2moderate3high);table(apu)
lexis_df3$low2moderate3high<-factor(apu,levels=1:3,labels =c("low","mid", "high"))


# Create a subgroup where low2moderate3high is either "low" or "mid"
low_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("low"), ]


mid_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("mid"), ]

high_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("high"), ]

tmpAF_FLULOW.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Unique_Drugs_ATC_Category,
                             data=low_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_low <- glm(cbind(lex.Xst=="IS",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia++age.c+sex+syopa+tulot3luokkaa+RenalVT+per.c,
                      data=low_subset,family="poisreg")


tmpAF_FLUMmoderate.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Unique_Drugs_ATC_Category,
                                   data=mid_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_moderate <- glm(cbind(lex.Xst=="IS",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                           data=mid_subset,family="poisreg")

tmpAF_FLUhigh.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Unique_Drugs_ATC_Category,
                              data=high_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_high <- glm(cbind(lex.Xst=="IS",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+stroke+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                       data=high_subset,family="poisreg")



sink(file="AdjustedAnalysislowmoderatehighUnique_Drugs_ATC_Category.html")

knitr::kable(round(ci.exp(tmp.AF_FLU_low),2),caption="low risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_moderate),2),caption="moderate risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_high),2),caption="high risk adjusted IRRs",format="html")

sink()

sink(file="RateTable_lowANDmoderateANDhighriskUnique_Drugs_ATC_Category.html")
knitr::kable(tmpAF_FLULOW.rt,caption="LOW Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUMmoderate.rt,caption="MODERATE Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUhigh.rt,caption="HIGH Rate table (1/100 person years)",format="html",digits=2)
sink()

#====================================================
# FULL PIPELINE: SPLINE IRR PLOT (FINAL)
#====================================================

library(splines)
library(ggplot2)
library(dplyr)
library(showtext)

# Font
font_add_google("Rosario", family = "rosario")
showtext_auto()

#----------------------------------------------------
# MODEL (Poisson with spline)
#----------------------------------------------------

tmp.m.spline <- glm(
  cbind(lex.Xst=="IS", lex.dur) ~ 
    ns(Unique_Drugs_ATC, df = 3) +
    Diabetes + per.c + Relevel(lex.Cst,list(1,2,3,4)) +
    sex + age.c + tulot3luokkaa + Hypertension + Hyperlipidemia +
    psykiatria + ALKO + dementia + syopa + anyvasc +
    Bleeding + stroke + LiverVT + RenalVT + HF,
  data = lexis_df3,
  family = "poisreg"
)

#----------------------------------------------------
# PREDICTION DATA
#----------------------------------------------------

newdata <- data.frame(
  Unique_Drugs_ATC = seq(
    min(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    max(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    length.out = 100
  )
)

newdata <- newdata %>%
  mutate(
    Diabetes = "no",
    per.c = levels(lexis_df3$per.c)[1],
    lex.Cst = levels(lexis_df3$lex.Cst)[1],
    sex = "male",
    age.c = levels(lexis_df3$age.c)[1],
    tulot3luokkaa = "low",
    Hypertension = "no",
    Hyperlipidemia = "no",
    psykiatria = "no",
    ALKO = "no",
    dementia = "no",
    syopa = "no",
    anyvasc = "no",
    Bleeding = "no",
    stroke = "no",
    LiverVT = "no",
    RenalVT = "no",
    HF = "no"
  )

newdata <- newdata %>%
  mutate(
    per.c = factor(per.c, levels = levels(lexis_df3$per.c)),
    lex.Cst = factor(lex.Cst, levels = levels(lexis_df3$lex.Cst)),
    age.c = factor(age.c, levels = levels(lexis_df3$age.c))
  )

#----------------------------------------------------
# PREDICT IRR (reference = 0 drugs)
#----------------------------------------------------

pred <- predict(tmp.m.spline, newdata = newdata, se.fit = TRUE)

ref_data <- newdata
ref_data$Unique_Drugs_ATC <- 0
ref_pred <- predict(tmp.m.spline, newdata = ref_data, se.fit = TRUE)

newdata$IRR <- exp(pred$fit - ref_pred$fit)
newdata$lower <- exp((pred$fit - 1.96 * pred$se.fit) - ref_pred$fit)
newdata$upper <- exp((pred$fit + 1.96 * pred$se.fit) - ref_pred$fit)

pIS <- ggplot(newdata, aes(x = Unique_Drugs_ATC, y = IRR)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#588692", alpha = 0.2) +
  
  geom_line(linewidth = 1, color = "#2C3E3F") +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             linewidth = 0.8,
             color = "black") +
  
  coord_cartesian(xlim = c(0, 20), ylim = c(0.25, 3)) +
  
  labs(
    x = "Number of drugs",
    y = "Adjusted IRR",
  ) +
  
  theme_classic(base_size = 14, base_family = "rosario") +
  
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, color = "black"),
    
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = grid::unit(0.2, "cm"),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0)
  ) +
  
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0))

# Show
pIS

#------------------------------
# SAVE (same as your KM workflow)
#------------------------------

ggsave("IRR_spline_plotppIS.pdf",
       plot = pIS,
       width = 7,
       height = 5.5,
       device = cairo_pdf,
       bg = "white")

pdf_image <- magick::image_read_pdf("IRR_spline_plotppIS.pdf", density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotppIS.png",
                    format = "png",
                    density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotpIS.tiff",
                    format = "tiff",
                    density = 300,
                    compression = "LZW")




####ICH####


library(haven)

cohort_df <- data.frame(read_sav(main_cohort_path))
head(cohort_df)



# Ensure date columns are in Date format
cohort_df$FirsteverICHAfterCohort <- as.Date(cohort_df$FirsteverICHAfterCohort)
cohort_df$RecurrentICHAfterCohort <- as.Date(cohort_df$RecurrentICHAfterCohort)
cohort_df$kuolpvmSPSSdate <- as.Date(cohort_df$kuolpvmSPSSdate)

# Step 1: ICHfirstdayaftercohortfirstorrecurrent
cohort_df$ICHfirstdayaftercohortfirstorrecurrent <- pmin(
  cohort_df$FirsteverICHAfterCohort,
  cohort_df$RecurrentICHAfterCohort,
  na.rm = TRUE
)

# Step 2: ENDdateICH_loppuICHordeath (include fixed 2018-12-31)
fixed_end_date <- as.Date("2018-12-31")

cohort_df$ENDdateICH_loppuICHordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$ICHfirstdayaftercohortfirstorrecurrent,
  fixed_end_date,
  na.rm = TRUE
)

cohort_df$ICHallaftercohort <- ifelse(
  !is.na(cohort_df$ICHfirstdayaftercohortfirstorrecurrent),
  1,
  0
)

#-----------------------------
# Load medication dataset
#-----------------------------
medication_df <- read_sav(medication_path)

# Quick preview
head(medication_df)

#-----------------------------
# Prepare medication data
#-----------------------------

# Ensure ATC codes are strings and create ATC level-4 variable
medication_df <- medication_df %>%
  dplyr::mutate(
    atc = as.character(atc),
    ATC_4 = substr(atc, 1, 4)
  )

# Keep only relevant columns and remove duplicates (unique ATC per SID)
medication_unique_df <- data %>%
  dplyr::select(SID, atc, ATC_4) %>%
  dplyr::distinct()

# Summarize number of unique drugs per patient
medication_summary_df <- medication_unique_df %>%
  dplyr::group_by(SID) %>%
  dplyr::summarise(
    Unique_Drugs_ATC = dplyr::n_distinct(atc),
    Unique_Drugs_ATC4 = dplyr::n_distinct(ATC_4),
    .groups = "drop"
  )

#-----------------------------
# Merge summaries into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(medication_summary_df, by = "SID") %>%
  dplyr::mutate(
    Unique_Drugs_ATC = ifelse(is.na(Unique_Drugs_ATC), 0, Unique_Drugs_ATC),
    Unique_Drugs_ATC4 = ifelse(is.na(Unique_Drugs_ATC4), 0, Unique_Drugs_ATC4)
  )

#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC, Unique_Drugs_ATC4) %>%
  head()


library(dplyr)

# Categorize Unique_Drugs_ATC into 4 groups
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = dplyr::case_when(
      Unique_Drugs_ATC >= 0  & Unique_Drugs_ATC <= 2  ~ "0–2 drugs",
      Unique_Drugs_ATC >= 3  & Unique_Drugs_ATC <= 4  ~ "3–4 drugs",
      Unique_Drugs_ATC >= 5  & Unique_Drugs_ATC <= 9  ~ "5–9 drugs",
      Unique_Drugs_ATC >= 10                          ~ "10+ drugs",
      TRUE                                            ~ NA_character_   # fallback
    )
  )

# Optional: make it an ordered factor
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = factor(
      Unique_Drugs_ATC_Category,
      levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
      ordered = TRUE
    )
  )

# Quick check
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_Category,
  levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)

with(cohort_df,range(cal.yr(ENDdateICH_loppuICHordeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","ICH")))
exit_status_factor<-with(cohort_df,factor(ifelse(ICHallaftercohort==1,"ICH","EOF"),levels=c("SOF","EOF","ICH")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(ENDdateICH_loppuICHordeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)

library(cohorttools)

boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))


library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="ICH")~Unique_Drugs_ATC_Category,
                    data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTableUnique_Drugs_ATC_CategoryICH.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="ICH",lex.dur) ~ ALKO+ Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiomICH <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiomICH <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiomICH <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)

sink(file="IRRTulosUnique_Drugs_ATC_CategoryICH.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()






apu<-unclass(lexis_df3$low2moderate3high);table(apu)
lexis_df3$low2moderate3high<-factor(apu,levels=1:3,labels =c("low","mid", "high"))


# Create a subgroup where low2moderate3high is either "low" or "mid"
low_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("low"), ]


mid_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("mid"), ]

high_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("high"), ]

tmpAF_FLULOW.rt<-mkratetable(Surv(lex.dur,lex.Xst=="ICH")~Unique_Drugs_ATC_Category,
                             data=low_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_low <- glm(cbind(lex.Xst=="ICH",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia++age.c+sex+syopa+tulot3luokkaa+RenalVT+per.c,
                      data=low_subset,family="poisreg")


tmpAF_FLUMmoderate.rt<-mkratetable(Surv(lex.dur,lex.Xst=="ICH")~Unique_Drugs_ATC_Category,
                                   data=mid_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_moderate <- glm(cbind(lex.Xst=="ICH",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                           data=mid_subset,family="poisreg")

tmpAF_FLUhigh.rt<-mkratetable(Surv(lex.dur,lex.Xst=="ICH")~Unique_Drugs_ATC_Category,
                              data=high_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_high <- glm(cbind(lex.Xst=="ICH",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+stroke+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                       data=high_subset,family="poisreg")



sink(file="AdjustedAnalysislowmoderatehighUnique_Drugs_ATC_CategoryICH.html")

knitr::kable(round(ci.exp(tmp.AF_FLU_low),2),caption="low risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_moderate),2),caption="moderate risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_high),2),caption="high risk adjusted IRRs",format="html")

sink()

sink(file="RateTable_lowANDmoderateANDhighriskUnique_Drugs_ATC_CategoryICH.html")
knitr::kable(tmpAF_FLULOW.rt,caption="LOW Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUMmoderate.rt,caption="MODERATE Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUhigh.rt,caption="HIGH Rate table (1/100 person years)",format="html",digits=2)
sink()


#----------------------------------------------------
# MODEL (Poisson with spline)
#----------------------------------------------------

tmp.m.spline <- glm(
  cbind(lex.Xst=="ICH", lex.dur) ~ 
    ns(Unique_Drugs_ATC, df = 3) +
    Diabetes + per.c + Relevel(lex.Cst,list(1,2,3,4)) +
    sex + age.c + tulot3luokkaa + Hypertension + Hyperlipidemia +
    psykiatria + ALKO + dementia + syopa + anyvasc +
    Bleeding + stroke + LiverVT + RenalVT + HF,
  data = lexis_df3,
  family = "poisreg"
)

#----------------------------------------------------
# PREDICTION DATA
#----------------------------------------------------

newdata <- data.frame(
  Unique_Drugs_ATC = seq(
    min(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    max(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    length.out = 100
  )
)

newdata <- newdata %>%
  mutate(
    Diabetes = "no",
    per.c = levels(lexis_df3$per.c)[1],
    lex.Cst = levels(lexis_df3$lex.Cst)[1],
    sex = "male",
    age.c = levels(lexis_df3$age.c)[1],
    tulot3luokkaa = "low",
    Hypertension = "no",
    Hyperlipidemia = "no",
    psykiatria = "no",
    ALKO = "no",
    dementia = "no",
    syopa = "no",
    anyvasc = "no",
    Bleeding = "no",
    stroke = "no",
    LiverVT = "no",
    RenalVT = "no",
    HF = "no"
  )

newdata <- newdata %>%
  mutate(
    per.c = factor(per.c, levels = levels(lexis_df3$per.c)),
    lex.Cst = factor(lex.Cst, levels = levels(lexis_df3$lex.Cst)),
    age.c = factor(age.c, levels = levels(lexis_df3$age.c))
  )

#----------------------------------------------------
# PREDICT IRR (reference = 0 drugs)
#----------------------------------------------------

pred <- predict(tmp.m.spline, newdata = newdata, se.fit = TRUE)

ref_data <- newdata
ref_data$Unique_Drugs_ATC <- 0
ref_pred <- predict(tmp.m.spline, newdata = ref_data, se.fit = TRUE)

newdata$IRR <- exp(pred$fit - ref_pred$fit)
newdata$lower <- exp((pred$fit - 1.96 * pred$se.fit) - ref_pred$fit)
newdata$upper <- exp((pred$fit + 1.96 * pred$se.fit) - ref_pred$fit)

pICH <- ggplot(newdata, aes(x = Unique_Drugs_ATC, y = IRR)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#588692", alpha = 0.2) +
  
  geom_line(linewidth = 1, color = "#2C3E3F") +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             linewidth = 0.8,
             color = "black") +
  
  coord_cartesian(xlim = c(0, 20), ylim = c(0.25, 3)) +
  
  labs(
    x = "Number of drugs",
    y = "Adjusted IRR",
  ) +
  
  theme_classic(base_size = 14, base_family = "rosario") +
  
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, color = "black"),
    
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = grid::unit(0.2, "cm"),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0)
  ) +
  
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0))

# Show
pICH

#------------------------------
# SAVE (same as your KM workflow)
#------------------------------

ggsave("IRR_spline_plotpICH.pdf",
       plot = pICH,
       width = 7,
       height = 5.5,
       device = cairo_pdf,
       bg = "white")

pdf_image <- magick::image_read_pdf("IRR_spline_plotpICH.pdf", density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotpICH.png",
                    format = "png",
                    density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotpICH.tiff",
                    format = "tiff",
                    density = 300,
                    compression = "LZW")



######mortality#####
library(haven)

cohort_df <- data.frame(read_sav(main_cohort_path))
head(cohort_df)


#-----------------------------
# Load medication dataset
#-----------------------------
medication_df <- read_sav(medication_path)

# Quick preview
head(medication_df)

#-----------------------------
# Prepare medication data
#-----------------------------

# Ensure ATC codes are strings and create ATC level-4 variable
medication_df <- medication_df %>%
  dplyr::mutate(
    atc = as.character(atc),
    ATC_4 = substr(atc, 1, 4)
  )

# Keep only relevant columns and remove duplicates (unique ATC per SID)
medication_unique_df <- data %>%
  dplyr::select(SID, atc, ATC_4) %>%
  dplyr::distinct()

# Summarize number of unique drugs per patient
medication_summary_df <- medication_unique_df %>%
  dplyr::group_by(SID) %>%
  dplyr::summarise(
    Unique_Drugs_ATC = dplyr::n_distinct(atc),
    Unique_Drugs_ATC4 = dplyr::n_distinct(ATC_4),
    .groups = "drop"
  )

#-----------------------------
# Merge summaries into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(medication_summary_df, by = "SID") %>%
  dplyr::mutate(
    Unique_Drugs_ATC = ifelse(is.na(Unique_Drugs_ATC), 0, Unique_Drugs_ATC),
    Unique_Drugs_ATC4 = ifelse(is.na(Unique_Drugs_ATC4), 0, Unique_Drugs_ATC4)
  )

#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC, Unique_Drugs_ATC4) %>%
  head()


library(dplyr)

# Categorize Unique_Drugs_ATC into 4 groups
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = dplyr::case_when(
      Unique_Drugs_ATC >= 0  & Unique_Drugs_ATC <= 2  ~ "0–2 drugs",
      Unique_Drugs_ATC >= 3  & Unique_Drugs_ATC <= 4  ~ "3–4 drugs",
      Unique_Drugs_ATC >= 5  & Unique_Drugs_ATC <= 9  ~ "5–9 drugs",
      Unique_Drugs_ATC >= 10                          ~ "10+ drugs",
      TRUE                                            ~ NA_character_   # fallback
    )
  )

# Optional: make it an ordered factor
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = factor(
      Unique_Drugs_ATC_Category,
      levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
      ordered = TRUE
    )
  )

# Quick check
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_Category,
  levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)

with(cohort_df,range(cal.yr(DateLoppuOrDeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","Dead")))
exit_status_factor<-with(cohort_df,factor(ifelse(Kuollut==1,"Dead","EOF"),levels=c("SOF","EOF","Dead")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateLoppuOrDeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)

library(cohorttools)

boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))



apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))



library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="Dead")~per.c+age.c+Relevel(lex.Cst,list(1,2,3,4))+
                      sex+Hyperlipidemia+Hypertension+Diabetes+psykiatria+dementia+syopa+anyvasc+DM1onlyins2insoral3oral4noMed+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
                    data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTableUnique_Drugs_ATC_Categorymortality.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="Dead",lex.dur) ~ ALKO+Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiomDead <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiomDead <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiomDead <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)

sink(file="IRRTulosUnique_Drugs_ATC_Categorymortality.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()




#####


apu<-unclass(lexis_df3$low2moderate3high);table(apu)
lexis_df3$low2moderate3high<-factor(apu,levels=1:3,labels =c("low","mid", "high"))


# Create a subgroup where low2moderate3high is either "low" or "mid"
low_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("low"), ]


mid_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("mid"), ]

high_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("high"), ]

tmpAF_FLULOW.rt<-mkratetable(Surv(lex.dur,lex.Xst=="Dead")~Unique_Drugs_ATC_Category,
                             data=low_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_low <- glm(cbind(lex.Xst=="Dead",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia++age.c+sex+syopa+tulot3luokkaa+RenalVT+per.c,
                      data=low_subset,family="poisreg")


tmpAF_FLUMmoderate.rt<-mkratetable(Surv(lex.dur,lex.Xst=="Dead")~Unique_Drugs_ATC_Category,
                                   data=mid_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_moderate <- glm(cbind(lex.Xst=="Dead",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                           data=mid_subset,family="poisreg")

tmpAF_FLUhigh.rt<-mkratetable(Surv(lex.dur,lex.Xst=="Dead")~Unique_Drugs_ATC_Category,
                              data=high_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_high <- glm(cbind(lex.Xst=="Dead",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+stroke+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                       data=high_subset,family="poisreg")



sink(file="AdjustedAnalysislowmoderatehighUnique_Drugs_ATC_CategoryDead.html")

knitr::kable(round(ci.exp(tmp.AF_FLU_low),2),caption="low risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_moderate),2),caption="moderate risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_high),2),caption="high risk adjusted IRRs",format="html")

sink()

sink(file="RateTable_lowANDmoderateANDhighriskUnique_Drugs_ATC_CategoryDead.html")
knitr::kable(tmpAF_FLULOW.rt,caption="LOW Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUMmoderate.rt,caption="MODERATE Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUhigh.rt,caption="HIGH Rate table (1/100 person years)",format="html",digits=2)
sink()

#====================================================
# FULL PIPELINE: SPLINE IRR PLOT (FINAL)
#====================================================

library(splines)
library(ggplot2)
library(dplyr)
library(showtext)

# Font
font_add_google("Rosario", family = "rosario")
showtext_auto()

#----------------------------------------------------
# MODEL (Poisson with spline)
#----------------------------------------------------

tmp.m.spline <- glm(
  cbind(lex.Xst=="Dead", lex.dur) ~ 
    ns(Unique_Drugs_ATC, df = 3) +
    Diabetes + per.c + Relevel(lex.Cst,list(1,2,3,4)) +
    sex + age.c + tulot3luokkaa + Hypertension + Hyperlipidemia +
    psykiatria + ALKO + dementia + syopa + anyvasc +
    Bleeding + stroke + LiverVT + RenalVT + HF,
  data = lexis_df3,
  family = "poisreg"
)

#----------------------------------------------------
# PREDICTION DATA
#----------------------------------------------------

newdata <- data.frame(
  Unique_Drugs_ATC = seq(
    min(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    max(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    length.out = 100
  )
)

newdata <- newdata %>%
  mutate(
    Diabetes = "no",
    per.c = levels(lexis_df3$per.c)[1],
    lex.Cst = levels(lexis_df3$lex.Cst)[1],
    sex = "male",
    age.c = levels(lexis_df3$age.c)[1],
    tulot3luokkaa = "low",
    Hypertension = "no",
    Hyperlipidemia = "no",
    psykiatria = "no",
    ALKO = "no",
    dementia = "no",
    syopa = "no",
    anyvasc = "no",
    Bleeding = "no",
    stroke = "no",
    LiverVT = "no",
    RenalVT = "no",
    HF = "no"
  )

newdata <- newdata %>%
  mutate(
    per.c = factor(per.c, levels = levels(lexis_df3$per.c)),
    lex.Cst = factor(lex.Cst, levels = levels(lexis_df3$lex.Cst)),
    age.c = factor(age.c, levels = levels(lexis_df3$age.c))
  )

#----------------------------------------------------
# PREDICT IRR (reference = 0 drugs)
#----------------------------------------------------

pred <- predict(tmp.m.spline, newdata = newdata, se.fit = TRUE)

ref_data <- newdata
ref_data$Unique_Drugs_ATC <- 0
ref_pred <- predict(tmp.m.spline, newdata = ref_data, se.fit = TRUE)

newdata$IRR <- exp(pred$fit - ref_pred$fit)
newdata$lower <- exp((pred$fit - 1.96 * pred$se.fit) - ref_pred$fit)
newdata$upper <- exp((pred$fit + 1.96 * pred$se.fit) - ref_pred$fit)

pDead <- ggplot(newdata, aes(x = Unique_Drugs_ATC, y = IRR)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#588692", alpha = 0.2) +
  
  geom_line(linewidth = 1, color = "#2C3E3F") +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             linewidth = 0.8,
             color = "black") +
  
  coord_cartesian(xlim = c(0, 20), ylim = c(0.25, 3)) +
  
  labs(
    x = "Number of drugs",
    y = "Adjusted IRR",
  ) +
  
  theme_classic(base_size = 14, base_family = "rosario") +
  
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, color = "black"),
    
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = grid::unit(0.2, "cm"),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0)
  ) +
  
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0))

# Show
pDead

#------------------------------
# SAVE (same as your KM workflow)
#------------------------------

ggsave("IRR_spline_plotppDead.pdf",
       plot = pDead,
       width = 7,
       height = 5.5,
       device = cairo_pdf,
       bg = "white")

pdf_image <- magick::image_read_pdf("IRR_spline_plotppDead.pdf", density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotppDead.png",
                    format = "png",
                    density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotpDead.tiff",
                    format = "tiff",
                    density = 300,
                    compression = "LZW")




######any bleeding #####

cohort_df <- data.frame(read_sav(main_cohort_path))



# Ensure date columns are in Date format
cohort_df$CombinedFirstEverRecurrentAfterCohort <- as.Date(cohort_df$CombinedFirstEverRecurrentAfterCohort)

cohort_df$kuolpvmSPSSdate <- as.Date(cohort_df$kuolpvmSPSSdate)



# Step 2: ENDdateICH_loppuICHordeath (include fixed 2018-12-31)
fixed_end_date <- as.Date("2018-12-31")

cohort_df$ENDdateAnyBleed_loppuAnyBleedordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$CombinedFirstEverRecurrentAfterCohort,
  fixed_end_date,
  na.rm = TRUE
)

cohort_df$Anybleedallaftercohort <- ifelse(
  !is.na(cohort_df$CombinedFirstEverRecurrentAfterCohort),
  1,
  0
)




#-----------------------------
# Load medication dataset
#-----------------------------
medication_df <- read_sav(medication_path)

# Quick preview
head(medication_df)

#-----------------------------
# Prepare medication data
#-----------------------------

# Ensure ATC codes are strings and create ATC level-4 variable
medication_df <- medication_df %>%
  dplyr::mutate(
    atc = as.character(atc),
    ATC_4 = substr(atc, 1, 4)
  )

# Keep only relevant columns and remove duplicates (unique ATC per SID)
medication_unique_df <- data %>%
  dplyr::select(SID, atc, ATC_4) %>%
  dplyr::distinct()

# Summarize number of unique drugs per patient
medication_summary_df <- medication_unique_df %>%
  dplyr::group_by(SID) %>%
  dplyr::summarise(
    Unique_Drugs_ATC = dplyr::n_distinct(atc),
    Unique_Drugs_ATC4 = dplyr::n_distinct(ATC_4),
    .groups = "drop"
  )

#-----------------------------
# Merge summaries into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(medication_summary_df, by = "SID") %>%
  dplyr::mutate(
    Unique_Drugs_ATC = ifelse(is.na(Unique_Drugs_ATC), 0, Unique_Drugs_ATC),
    Unique_Drugs_ATC4 = ifelse(is.na(Unique_Drugs_ATC4), 0, Unique_Drugs_ATC4)
  )

#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC, Unique_Drugs_ATC4) %>%
  head()


library(dplyr)

# Categorize Unique_Drugs_ATC into 4 groups
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = dplyr::case_when(
      Unique_Drugs_ATC >= 0  & Unique_Drugs_ATC <= 2  ~ "0–2 drugs",
      Unique_Drugs_ATC >= 3  & Unique_Drugs_ATC <= 4  ~ "3–4 drugs",
      Unique_Drugs_ATC >= 5  & Unique_Drugs_ATC <= 9  ~ "5–9 drugs",
      Unique_Drugs_ATC >= 10                          ~ "10+ drugs",
      TRUE                                            ~ NA_character_   # fallback
    )
  )

# Optional: make it an ordered factor
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = factor(
      Unique_Drugs_ATC_Category,
      levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
      ordered = TRUE
    )
  )

# Quick check
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_Category,
  levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)

with(cohort_df,range(cal.yr(ENDdateAnyBleed_loppuAnyBleedordeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","AnyBleed")))
exit_status_factor<-with(cohort_df,factor(ifelse(Anybleedallaftercohort==1,"AnyBleed","EOF"),levels=c("SOF","EOF","AnyBleed")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(ENDdateAnyBleed_loppuAnyBleedordeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)

library(cohorttools)

boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))



apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))



library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="AnyBleed")~per.c+age.c+Relevel(lex.Cst,list(1,2,3,4))+
                      sex+Hyperlipidemia+Hypertension+Diabetes+psykiatria+dementia+syopa+anyvasc+DM1onlyins2insoral3oral4noMed+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
                    data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTableUnique_Drugs_ATC_CategoryAnyBleed.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="AnyBleed",lex.dur) ~ ALKO+Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiomAnyBleed <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiomAnyBleed <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiomAnyBleed <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)

sink(file="IRRTulosUnique_Drugs_ATC_CategoryAnyBleed.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()




apu<-unclass(lexis_df3$low2moderate3high);table(apu)
lexis_df3$low2moderate3high<-factor(apu,levels=1:3,labels =c("low","mid", "high"))


# Create a subgroup where low2moderate3high is either "low" or "mid"
low_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("low"), ]


mid_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("mid"), ]

high_subset <- lexis_df3[lexis_df3$low2moderate3high %in% c("high"), ]

tmpAF_FLULOW.rt<-mkratetable(Surv(lex.dur,lex.Xst=="AnyBleed")~Unique_Drugs_ATC_Category,
                             data=low_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_low <- glm(cbind(lex.Xst=="AnyBleed",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia++age.c+sex+syopa+tulot3luokkaa+RenalVT+per.c,
                      data=low_subset,family="poisreg")


tmpAF_FLUMmoderate.rt<-mkratetable(Surv(lex.dur,lex.Xst=="AnyBleed")~Unique_Drugs_ATC_Category,
                                   data=mid_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_moderate <- glm(cbind(lex.Xst=="AnyBleed",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                           data=mid_subset,family="poisreg")

tmpAF_FLUhigh.rt<-mkratetable(Surv(lex.dur,lex.Xst=="AnyBleed")~Unique_Drugs_ATC_Category,
                              data=high_subset,scale=100,add.RR = TRUE)

tmp.AF_FLU_high <- glm(cbind(lex.Xst=="AnyBleed",lex.dur) ~ ALKO+Unique_Drugs_ATC_Category+Relevel(lex.Cst,list(1,2,3,4))+Hyperlipidemia+Diabetes++age.c+sex+syopa+stroke+tulot3luokkaa+Hypertension+anyvasc+RenalVT+HF+per.c,
                       data=high_subset,family="poisreg")



sink(file="AdjustedAnalysislowmoderatehighUnique_Drugs_ATC_CategoryAnyBleed.html")

knitr::kable(round(ci.exp(tmp.AF_FLU_low),2),caption="low risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_moderate),2),caption="moderate risk adjusted IRRs",format="html")
knitr::kable(round(ci.exp(tmp.AF_FLU_high),2),caption="high risk adjusted IRRs",format="html")

sink()

sink(file="RateTable_lowANDmoderateANDhighriskUnique_Drugs_ATC_CategoryAnyBleed.html")
knitr::kable(tmpAF_FLULOW.rt,caption="LOW Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUMmoderate.rt,caption="MODERATE Rate table (1/100 person years)",format="html",digits=2)
knitr::kable(tmpAF_FLUhigh.rt,caption="HIGH Rate table (1/100 person years)",format="html",digits=2)
sink()


#====================================================
# FULL PIPELINE: SPLINE IRR PLOT (FINAL)
#====================================================

library(splines)
library(ggplot2)
library(dplyr)
library(showtext)

# Font
font_add_google("Rosario", family = "rosario")
showtext_auto()

#----------------------------------------------------
# MODEL (Poisson with spline)
#----------------------------------------------------

tmp.m.spline <- glm(
  cbind(lex.Xst=="AnyBleed", lex.dur) ~ 
    ns(Unique_Drugs_ATC, df = 3) +
    Diabetes + per.c + Relevel(lex.Cst,list(1,2,3,4)) +
    sex + age.c + tulot3luokkaa + Hypertension + Hyperlipidemia +
    psykiatria + ALKO + dementia + syopa + anyvasc +
    Bleeding + stroke + LiverVT + RenalVT + HF,
  data = lexis_df3,
  family = "poisreg"
)

#----------------------------------------------------
# PREDICTION DATA
#----------------------------------------------------

newdata <- data.frame(
  Unique_Drugs_ATC = seq(
    min(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    max(lexis_df3$Unique_Drugs_ATC, na.rm = TRUE),
    length.out = 100
  )
)

newdata <- newdata %>%
  mutate(
    Diabetes = "no",
    per.c = levels(lexis_df3$per.c)[1],
    lex.Cst = levels(lexis_df3$lex.Cst)[1],
    sex = "male",
    age.c = levels(lexis_df3$age.c)[1],
    tulot3luokkaa = "low",
    Hypertension = "no",
    Hyperlipidemia = "no",
    psykiatria = "no",
    ALKO = "no",
    dementia = "no",
    syopa = "no",
    anyvasc = "no",
    Bleeding = "no",
    stroke = "no",
    LiverVT = "no",
    RenalVT = "no",
    HF = "no"
  )

newdata <- newdata %>%
  mutate(
    per.c = factor(per.c, levels = levels(lexis_df3$per.c)),
    lex.Cst = factor(lex.Cst, levels = levels(lexis_df3$lex.Cst)),
    age.c = factor(age.c, levels = levels(lexis_df3$age.c))
  )

#----------------------------------------------------
# PREDICT IRR (reference = 0 drugs)
#----------------------------------------------------

pred <- predict(tmp.m.spline, newdata = newdata, se.fit = TRUE)

ref_data <- newdata
ref_data$Unique_Drugs_ATC <- 0
ref_pred <- predict(tmp.m.spline, newdata = ref_data, se.fit = TRUE)

newdata$IRR <- exp(pred$fit - ref_pred$fit)
newdata$lower <- exp((pred$fit - 1.96 * pred$se.fit) - ref_pred$fit)
newdata$upper <- exp((pred$fit + 1.96 * pred$se.fit) - ref_pred$fit)

pAnyBleed <- ggplot(newdata, aes(x = Unique_Drugs_ATC, y = IRR)) +
  
  geom_ribbon(aes(ymin = lower, ymax = upper),
              fill = "#588692", alpha = 0.2) +
  
  geom_line(linewidth = 1, color = "#2C3E3F") +
  
  geom_hline(yintercept = 1,
             linetype = "dashed",
             linewidth = 0.8,
             color = "black") +
  
  coord_cartesian(xlim = c(0, 20), ylim = c(0.25, 3)) +
  
  labs(
    x = "Number of drugs",
    y = "Adjusted IRR",
  ) +
  
  theme_classic(base_size = 14, base_family = "rosario") +
  
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 13, color = "black"),
    
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = grid::unit(0.2, "cm"),
    
    plot.title = element_text(size = 16, face = "bold", hjust = 0)
  ) +
  
  scale_x_continuous(breaks = seq(0, 20, by = 5)) +
  scale_y_continuous(expand = c(0, 0))

# Show
pAnyBleed

#------------------------------
# SAVE (same as your KM workflow)
#------------------------------

ggsave("IRR_spline_plotppAnyBleed.pdf",
       plot = pAnyBleed,
       width = 7,
       height = 5.5,
       device = cairo_pdf,
       bg = "white")

pdf_image <- magick::image_read_pdf("IRR_spline_plotppAnyBleed.pdf", density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotppAnyBleed.png",
                    format = "png",
                    density = 300)

magick::image_write(pdf_image,
                    path = "IRR_spline_plotpAnyBleed.tiff",
                    format = "tiff",
                    density = 300,
                    compression = "LZW")

# =====================================================
# FOLLOW-UP TIME SUMMARY (per outcome)
# =====================================================

calc_fu <- function(x) {
  c(
    mean = mean(x, na.rm = TRUE),
    sd   = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    n = length(x)
  )
}

fu_results <- rbind(
  Stroke = calc_fu(tmp_stroke$duration_stroke),
  Mortality = calc_fu(tmp_mortality$duration_mortality),
  ICH = calc_fu(tmp_ich$duration_ich),
  MajorBleed = calc_fu(tmp_bleed$duration_bleed)
)

fu_results <- round(fu_results, 2)

print(fu_results)

# =====================================================
# OVERALL FOLLOW-UP (ALL OUTCOMES)
# =====================================================

cohort_df$duration_overall <- as.numeric(
  difftime(cohort_df$DateLoppuOrDeath,
           cohort_df$CohortEntryDate,
           units = "days")
) / 365.25

tmp_overall <- cohort_df %>%
  dplyr::filter(!is.na(duration_overall) & duration_overall >= 0)

overall_fu <- c(
  mean = mean(tmp_overall$duration_overall, na.rm = TRUE),
  sd   = sd(tmp_overall$duration_overall, na.rm = TRUE),
  median = median(tmp_overall$duration_overall, na.rm = TRUE),
  n = sum(!is.na(tmp_overall$duration_overall))
)

round(overall_fu, 2)

####### plots ######

# --- Load libraries ---
library(haven)
library(dplyr)
library(lubridate)
library(survival)
library(survminer)
library(showtext)
library(sysfonts)
library(cowplot)
library(magick)
# =====================================================
# --- Load data and libraries ---
# =====================================================
library(haven)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(cowplot)
library(showtext)
library(magick)
library(sysfonts)

# --- Read dataset ---
cohort_df <- data.frame(read_sav(main_cohort_path))


#-----------------------------
# Load medication dataset
#-----------------------------
medication_df <- read_sav(medication_path)

# Quick preview
head(medication_df)

#-----------------------------
# Prepare medication data
#-----------------------------

# Ensure ATC codes are strings and create ATC level-4 variable
medication_df <- medication_df %>%
  dplyr::mutate(
    atc = as.character(atc),
    ATC_4 = substr(atc, 1, 4)
  )

# Keep only relevant columns and remove duplicates (unique ATC per SID)
medication_unique_df <- data %>%
  dplyr::select(SID, atc, ATC_4) %>%
  dplyr::distinct()

# Summarize number of unique drugs per patient
medication_summary_df <- medication_unique_df %>%
  dplyr::group_by(SID) %>%
  dplyr::summarise(
    Unique_Drugs_ATC = dplyr::n_distinct(atc),
    Unique_Drugs_ATC4 = dplyr::n_distinct(ATC_4),
    .groups = "drop"
  )

#-----------------------------
# Merge summaries into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(medication_summary_df, by = "SID") %>%
  dplyr::mutate(
    Unique_Drugs_ATC = ifelse(is.na(Unique_Drugs_ATC), 0, Unique_Drugs_ATC),
    Unique_Drugs_ATC4 = ifelse(is.na(Unique_Drugs_ATC4), 0, Unique_Drugs_ATC4)
  )

#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC, Unique_Drugs_ATC4) %>%
  head()


library(dplyr)

# Categorize Unique_Drugs_ATC into 4 groups
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = dplyr::case_when(
      Unique_Drugs_ATC >= 0  & Unique_Drugs_ATC <= 2  ~ "0–2 drugs",
      Unique_Drugs_ATC >= 3  & Unique_Drugs_ATC <= 4  ~ "3–4 drugs",
      Unique_Drugs_ATC >= 5  & Unique_Drugs_ATC <= 9  ~ "5–9 drugs",
      Unique_Drugs_ATC >= 10                          ~ "10+ drugs",
      TRUE                                            ~ NA_character_   # fallback
    )
  )

# Optional: make it an ordered factor
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = factor(
      Unique_Drugs_ATC_Category,
      levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
      ordered = TRUE
    )
  )

# Quick check
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)


# --- Exposure variable: Unique_Drugs_ATC_Category (already created earlier) ---
# If not created yet, you can include the categorization here:
cohort_df <- cohort_df %>%
  mutate(
    Unique_Drugs_ATC_Category = case_when(
      Unique_Drugs_ATC >= 0 & Unique_Drugs_ATC <= 2 ~ "0–2 drugs",
      Unique_Drugs_ATC >= 3 & Unique_Drugs_ATC <= 4 ~ "3–4 drugs",
      Unique_Drugs_ATC >= 5 & Unique_Drugs_ATC <= 9 ~ "5–9 drugs",
      Unique_Drugs_ATC >= 10 ~ "10+ drugs",
      TRUE ~ NA_character_
    ),
    Unique_Drugs_ATC_Category = factor(
      Unique_Drugs_ATC_Category,
      levels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
      ordered = TRUE
    )
  )

# --- Font setup ---
font_add_google("Rosario", "rosario")
showtext_auto()

# --- Fixed study end date ---
fixed_end_date <- as.Date("2018-12-31")

# --- Ensure all relevant date variables are in Date format ---
date_vars <- c("CohortEntryDate",
               "DateISorLoppuOrDeath",
               "DateLoppuOrDeath",
               "FirsteverICHAfterCohort",
               "RecurrentICHAfterCohort",
               "kuolpvmSPSSdate",
               "CombinedFirstEverRecurrentAfterCohort")

for (v in date_vars) {
  if (v %in% names(cohort_df)) {
    cohort_df[[v]] <- as.Date(cohort_df[[v]])
  }
}

# =====================================================
# 1. STROKE
# =====================================================
cohort_df$duration_stroke <- as.numeric(
  difftime(cohort_df$DateISorLoppuOrDeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25

tmp_stroke <- cohort_df %>% filter(!is.na(duration_stroke) & duration_stroke >= 0)

surv_stroke <- with(tmp_stroke, Surv(duration_stroke, ISaftercohortall == 1))
km_stroke <- survfit(surv_stroke ~ Unique_Drugs_ATC_Category, data = tmp_stroke)

# =====================================================
# 2. MORTALITY
# =====================================================
cohort_df$duration_mortality <- as.numeric(
  difftime(cohort_df$DateLoppuOrDeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25

tmp_mortality <- cohort_df %>% filter(!is.na(duration_mortality) & duration_mortality >= 0)

surv_mortality <- with(tmp_mortality, Surv(duration_mortality, Kuollut == 1))
km_mortality <- survfit(surv_mortality ~ Unique_Drugs_ATC_Category, data = tmp_mortality)

# =====================================================
# 3. INTRACRANIAL HEMORRHAGE (ICH)
# =====================================================
cohort_df$ICHfirstdayaftercohortfirstorrecurrent <- pmin(
  cohort_df$FirsteverICHAfterCohort,
  cohort_df$RecurrentICHAfterCohort,
  na.rm = TRUE
)

cohort_df$ENDdateICH_loppuICHordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$ICHfirstdayaftercohortfirstorrecurrent,
  fixed_end_date,
  na.rm = TRUE
)

cohort_df$ICHallaftercohort <- ifelse(
  !is.na(cohort_df$ICHfirstdayaftercohortfirstorrecurrent), 1, 0
)

cohort_df$duration_ich <- as.numeric(
  difftime(cohort_df$ENDdateICH_loppuICHordeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25

tmp_ich <- cohort_df %>% filter(!is.na(duration_ich) & duration_ich >= 0)

surv_ich <- with(tmp_ich, Surv(duration_ich, ICHallaftercohort == 1))
km_ich <- survfit(surv_ich ~ Unique_Drugs_ATC_Category, data = tmp_ich)

# =====================================================
# 4. MAJOR BLEEDING
# =====================================================
cohort_df$ENDdateAnyBleed_loppuAnyBleedordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$CombinedFirstEverRecurrentAfterCohort,
  fixed_end_date,
  na.rm = TRUE
)

cohort_df$Anybleedallaftercohort <- ifelse(
  !is.na(cohort_df$CombinedFirstEverRecurrentAfterCohort), 1, 0
)

cohort_df$duration_bleed <- as.numeric(
  difftime(cohort_df$ENDdateAnyBleed_loppuAnyBleedordeath,
           cohort_df$CohortEntryDate, units = "days")
) / 365.25

tmp_bleed <- cohort_df %>% filter(!is.na(duration_bleed) & duration_bleed >= 0)

surv_bleed <- with(tmp_bleed, Surv(duration_bleed, Anybleedallaftercohort == 1))
km_bleed <- survfit(surv_bleed ~ Unique_Drugs_ATC_Category, data = tmp_bleed)

# =====================================================
# Common Plot Settings
# =====================================================
# Use one color per drug category
custom_colors <- c("#588692", "#E6A960", "#B96D40", "#8C3F5D")

plot_fun <- function(fit, data, title_label, ylim_value) {
  p <- ggsurvplot(
    fit = fit,
    data = data,
    fun = "event",
    conf.int = TRUE,
    risk.table = FALSE,
    censor = FALSE,
    palette = c("#588692", "#E6A960", "#B96D40", "#8C3F5D"),  # optional custom palette
    pval = FALSE,
    xlab = "Follow-up time (years)",
    ylab = "Cumulative incidence",
    break.x.by = 1,
    legend.title = "Number of drugs",  # ✅ Clean title
    legend.labs = c("0–2", "3–4", "5–9", "10+"),  # ✅ No variable name prefix
    ggtheme = theme_classic(base_size = 14, base_family = "rosario") +
      theme(
        legend.title = element_text(size = 14, face = "bold"),
        legend.text = element_text(size = 13),
        axis.title = element_text(size = 16, face = "bold"),
        axis.text = element_text(size = 13),
        axis.ticks = element_line(linewidth = 0.8),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10))
      )
  )
  
  # Zoom to 5 years and adjust y-limits
  p$plot <- p$plot +
    coord_cartesian(xlim = c(0, 5), ylim = c(0, ylim_value)) +
    ggtitle(title_label)
  
  return(p)
}

# =====================================================
# Create plots A–D (no risk tables)
# =====================================================
plotA <- plot_fun(km_stroke, tmp_stroke, "A. Ischemic stroke", 0.15)
plotB <- plot_fun(km_mortality, tmp_mortality, "B. Mortality", 0.60)
plotC <- plot_fun(km_ich, tmp_ich, "C. Intracranial hemorrhage", 0.10)
plotD <- plot_fun(km_bleed, tmp_bleed, "D. Major bleeding", 0.30)

# =====================================================
# Combine into panel (no risk tables)
# =====================================================
panel_plot <- plot_grid(
  plotA$plot,
  plotB$plot,
  plotC$plot,
  plotD$plot,
  ncol = 2,
  labels = NULL
)

# ===================d==================================
# Save outputs
# =====================================================
ggsave("km_panel_Unique_Drugs_ATC_Category_no_risktable.pdf",
       plot = panel_plot, width = 16, height = 10)

pdf_image <- image_read_pdf("km_panel_Unique_Drugs_ATC_Category_no_risktable.pdf", density = 300)
image_write(pdf_image, path = "km_panel_Unique_Drugs_ATC_Category_no_risktable.png", format = "png", density = 300)
image_write(pdf_image, path = "km_panel_Unique_Drugs_ATC_Category_no_risktable.tiff", format = "tiff", density = 300, compression = "LZW")

print("✅ Kaplan–Meier panel (by Unique_Drugs_ATC_Category) saved as PDF, PNG, and TIFF.")

library(cowplot)
library(magick)

#-----------------------------------------------------
# ADD TITLES (MATCH KM STYLE)
#-----------------------------------------------------

pIS <- pIS +
  ggtitle("A. Ischemic stroke") +
  theme(
    plot.title = element_text(size = 14, face = "plain", hjust = 0)
  )

pDead <- pDead +
  ggtitle("B. Mortality") +
  theme(
    plot.title = element_text(size = 14, face = "plain", hjust = 0)
  )

pICH <- pICH +
  ggtitle("C. Intracranial hemorrhage") +
  theme(
    plot.title = element_text(size = 14, face = "plain", hjust = 0)
  )

pAnyBleed <- pAnyBleed +
  ggtitle("D. Major bleeding") +
  theme(
    plot.title = element_text(size = 14, face = "plain", hjust = 0)
  )

#-----------------------------------------------------
# COMBINE PANEL (2x2, NO EXTRA LABELS)
#-----------------------------------------------------

panel_plot <- plot_grid(
  pIS,
  pDead,
  pICH,
  pAnyBleed,
  ncol = 2
)

# Show
panel_plot

#-----------------------------------------------------
# SAVE (same workflow as KM panel)
#-----------------------------------------------------

ggsave("Spline_panel_IRR_2x2.pdf",
       plot = panel_plot,
       width = 16,
       height = 10,
       device = cairo_pdf,
       bg = "white")

pdf_image <- image_read_pdf("Spline_panel_IRR_2x2.pdf", density = 300)

image_write(pdf_image,
            path = "Spline_panel_IRR_2x2.png",
            format = "png",
            density = 300)

image_write(pdf_image,
            path = "Spline_panel_IRR_2x2.tiff",
            format = "tiff",
            density = 300,
            compression = "LZW")

print("✅ Spline panel (KM-style titles) saved as PDF, PNG, and TIFF.")


suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
  library(magick)
})

# Plain titles (no A/B/C/D prefixes), matching your chosen style
pIS <- pIS + ggtitle("A. Ischemic stroke")
pDead <- pDead + ggtitle("B. Mortality")
pICH <- pICH + ggtitle("C. Intracranial hemorrhage")
pAnyBleed <- pAnyBleed + ggtitle("D. Major bleeding")

# Match KM plain script typography and axis-label spacing
style_updates <- theme(
  axis.title = element_text(size = 17, face = "bold"),
  axis.text = element_text(size = 16, color = "black"),
  axis.ticks = element_line(linewidth = 0.8),
  axis.ticks.length = grid::unit(0.2, "cm"),
  axis.title.x = element_text(margin = margin(t = 18)),
  axis.title.y = element_text(margin = margin(r = 18)),
  plot.title = element_text(size = 17, face = "bold", hjust = 0)
)

pIS <- pIS + style_updates
pDead <- pDead + style_updates
pICH <- pICH + style_updates
pAnyBleed <- pAnyBleed + style_updates

# Same spacer-grid logic as KM panel script
null_plot <- ggplot() + theme_void()

panel_plot <- plot_grid(
  pIS,       null_plot, pDead,
  null_plot, null_plot, null_plot,
  pICH,      null_plot, pAnyBleed,
  ncol = 3,
  rel_widths = c(1, 0.07, 1),
  rel_heights = c(1, 0.07, 1),
  labels = NULL
)

# Save outputs
ggsave(
  "Spline_panel_IRR_2x2_o.pdf",
  plot = panel_plot,
  width = 16,
  height = 10,
  device = cairo_pdf,
  bg = "white"
)

pdf_image <- image_read_pdf("Spline_panel_IRR_2x2_o.pdf", density = 300)

image_write(
  pdf_image,
  path = "Spline_panel_IRR_2x2_o.png",
  format = "png",
  density = 300
)

image_write(
  pdf_image,
  path = "Spline_panel_IRR_2x2_o.tiff",
  format = "tiff",
  density = 300,
  compression = "LZW"
)

message("Wrote Spline_panel_IRR_2x2_o.pdf, .png, .tiff")

### Analysis with drug use after AF ###


######

# Load required packages
library(haven)
library(cohorttools)
library(dplyr)    # ensure dplyr functions are registered last

#-----------------------------
# Load main cohort dataset
#-----------------------------
cohort_df <- read_sav(main_cohort_path) %>%
  as.data.frame()

# Quick preview
head(cohort_df)
#-----------------------------
# Quick preview of cohort
#-----------------------------
head(cohort_df)

#-----------------------------
# Load post-AF polypharmacy category data
#-----------------------------
post_af_data <- read_sav(post_af_category_path) %>%
  as.data.frame() %>%
  dplyr::select(SID, Unique_Drugs_ATC_CategoryAFTER)

# Quick preview
head(post_af_data)

#-----------------------------
# Merge into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(post_af_data, by = "SID")

#-----------------------------
# Rename to match original variable name
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = Unique_Drugs_ATC_CategoryAFTER
  )

#-----------------------------
# Ensure correct factor structure
#-----------------------------
cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_CategoryAFTER,
  levels = c(1, 2, 3, 4),
  labels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)


#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC_Category) %>%
  head()


with(cohort_df,range(cal.yr(DateISorLoppuOrDeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","IS")))
exit_status_factor<-with(cohort_df,factor(ifelse(ISaftercohortall==1,"IS","EOF"),levels=c("SOF","EOF","IS")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateISorLoppuOrDeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)


boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))



library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="IS")~Unique_Drugs_ATC_Category, data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTablepolypharmacy_ISDRUGSAFTER.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="IS",lex.dur) ~ Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+
                sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+ALKO+
                dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiom <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiom <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiom <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)


sink(file="IRRTulosUnique_Drugs_ATC_Category_ISDRUGSAFTER.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()

####ICH####


library(haven)

cohort_df <- data.frame(read_sav(main_cohort_path))
head(cohort_df)



# Ensure date columns are in Date format
cohort_df$FirsteverICHAfterCohort <- as.Date(cohort_df$FirsteverICHAfterCohort)
cohort_df$RecurrentICHAfterCohort <- as.Date(cohort_df$RecurrentICHAfterCohort)
cohort_df$kuolpvmSPSSdate <- as.Date(cohort_df$kuolpvmSPSSdate)

# Step 1: ICHfirstdayaftercohortfirstorrecurrent
cohort_df$ICHfirstdayaftercohortfirstorrecurrent <- pmin(
  cohort_df$FirsteverICHAfterCohort,
  cohort_df$RecurrentICHAfterCohort,
  na.rm = TRUE
)

# Step 2: ENDdateICH_loppuICHordeath (include fixed 2018-12-31)
fixed_end_date <- as.Date("2018-12-31")

cohort_df$ENDdateICH_loppuICHordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$ICHfirstdayaftercohortfirstorrecurrent,
  fixed_end_date,
  na.rm = TRUE
)

cohort_df$ICHallaftercohort <- ifelse(
  !is.na(cohort_df$ICHfirstdayaftercohortfirstorrecurrent),
  1,
  0
)
#-----------------------------
# Quick preview of cohort
#-----------------------------
head(cohort_df)

#-----------------------------
# Load post-AF polypharmacy category data
#-----------------------------
post_af_data <- read_sav(post_af_category_path) %>%
  as.data.frame() %>%
  dplyr::select(SID, Unique_Drugs_ATC_CategoryAFTER)

# Quick preview
head(post_af_data)

#-----------------------------
# Merge into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(post_af_data, by = "SID")

#-----------------------------
# Rename to match original variable name
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = Unique_Drugs_ATC_CategoryAFTER
  )

#-----------------------------
# Ensure correct factor structure
#-----------------------------
cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_CategoryAFTER,
  levels = c(1, 2, 3, 4),
  labels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)


#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC_Category) %>%
  head()


with(cohort_df,range(cal.yr(ENDdateICH_loppuICHordeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","ICH")))
exit_status_factor<-with(cohort_df,factor(ifelse(ICHallaftercohort==1,"ICH","EOF"),levels=c("SOF","EOF","ICH")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(ENDdateICH_loppuICHordeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)

library(cohorttools)

boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))


library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="ICH")~Unique_Drugs_ATC_Category,
                    data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTableUnique_Drugs_ATC_CategoryICHDRUGSAFTER.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="ICH",lex.dur) ~ ALKO+ Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiomICH <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiomICH <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiomICH <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)

sink(file="IRRTulosUnique_Drugs_ATC_CategoryICHDRUGSAFTER.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()





######mortality#####
library(haven)

cohort_df <- data.frame(read_sav(main_cohort_path))
head(cohort_df)

#-----------------------------
# Quick preview of cohort
#-----------------------------
head(cohort_df)

#-----------------------------
# Load post-AF polypharmacy category data
#-----------------------------
post_af_data <- read_sav(post_af_category_path) %>%
  as.data.frame() %>%
  dplyr::select(SID, Unique_Drugs_ATC_CategoryAFTER)

# Quick preview
head(post_af_data)

#-----------------------------
# Merge into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(post_af_data, by = "SID")

#-----------------------------
# Rename to match original variable name
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = Unique_Drugs_ATC_CategoryAFTER
  )

#-----------------------------
# Ensure correct factor structure
#-----------------------------
cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_CategoryAFTER,
  levels = c(1, 2, 3, 4),
  labels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)


#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC_Category) %>%
  head()


with(cohort_df,range(cal.yr(DateLoppuOrDeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","Dead")))
exit_status_factor<-with(cohort_df,factor(ifelse(Kuollut==1,"Dead","EOF"),levels=c("SOF","EOF","Dead")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(DateLoppuOrDeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)

library(cohorttools)

boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))



apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))



library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="Dead")~per.c+age.c+Relevel(lex.Cst,list(1,2,3,4))+
                      sex+Hyperlipidemia+Hypertension+Diabetes+psykiatria+dementia+syopa+anyvasc+DM1onlyins2insoral3oral4noMed+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
                    data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTableUnique_Drugs_ATC_CategorymortalityDRUGSAFTER.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="Dead",lex.dur) ~ ALKO+Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiomDead <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiomDead <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiomDead <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)

sink(file="IRRTulosUnique_Drugs_ATC_CategorymortalityDRUGSAFTER.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()




######any bleeding #####

cohort_df <- data.frame(read_sav(main_cohort_path))



# Ensure date columns are in Date format
cohort_df$CombinedFirstEverRecurrentAfterCohort <- as.Date(cohort_df$CombinedFirstEverRecurrentAfterCohort)

cohort_df$kuolpvmSPSSdate <- as.Date(cohort_df$kuolpvmSPSSdate)



# Step 2: ENDdateICH_loppuICHordeath (include fixed 2018-12-31)
fixed_end_date <- as.Date("2018-12-31")

cohort_df$ENDdateAnyBleed_loppuAnyBleedordeath <- pmin(
  cohort_df$kuolpvmSPSSdate,
  cohort_df$CombinedFirstEverRecurrentAfterCohort,
  fixed_end_date,
  na.rm = TRUE
)

cohort_df$Anybleedallaftercohort <- ifelse(
  !is.na(cohort_df$CombinedFirstEverRecurrentAfterCohort),
  1,
  0
)


#-----------------------------
# Quick preview of cohort
#-----------------------------
head(cohort_df)

#-----------------------------
# Load post-AF polypharmacy category data
#-----------------------------
post_af_data <- read_sav(post_af_category_path) %>%
  as.data.frame() %>%
  dplyr::select(SID, Unique_Drugs_ATC_CategoryAFTER)

# Quick preview
head(post_af_data)

#-----------------------------
# Merge into main cohort
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::left_join(post_af_data, by = "SID")

#-----------------------------
# Rename to match original variable name
#-----------------------------
cohort_df <- cohort_df %>%
  dplyr::mutate(
    Unique_Drugs_ATC_Category = Unique_Drugs_ATC_CategoryAFTER
  )

#-----------------------------
# Ensure correct factor structure
#-----------------------------
cohort_df$Unique_Drugs_ATC_Category <- factor(
  cohort_df$Unique_Drugs_ATC_CategoryAFTER,
  levels = c(1, 2, 3, 4),
  labels = c("0–2 drugs", "3–4 drugs", "5–9 drugs", "10+ drugs"),
  ordered = FALSE
)


#-----------------------------
# Quick check of results
#-----------------------------
cohort_df %>%
  dplyr::count(Unique_Drugs_ATC_Category)

cohort_df %>%
  dplyr::select(SID, Unique_Drugs_ATC_Category) %>%
  head()


with(cohort_df,range(cal.yr(ENDdateAnyBleed_loppuAnyBleedordeath,format = "%Y-%m-%d") ))
with(cohort_df,range(cal.yr(kuolpvmSPSSdate,format = "%Y-%m-%d"),na.rm = TRUE))

entry_status_factor<-with(cohort_df,factor(rep("SOF",nrow(cohort_df)),levels=c("SOF","EOF","AnyBleed")))
exit_status_factor<-with(cohort_df,factor(ifelse(Anybleedallaftercohort==1,"AnyBleed","EOF"),levels=c("SOF","EOF","AnyBleed")))

table(entry_status_factor,useNA = "always")
table(exit_status_factor,useNA = "always")

lexis_df<-Lexis(entry=list(age=Age,
                         fu=0,
                         per=cal.yr(CohortEntryDate,format = "%Y-%m-%d")),
              duration = cal.yr(ENDdateAnyBleed_loppuAnyBleedordeath,format = "%Y-%m-%d")-
                cal.yr(CohortEntryDate,format = "%Y-%m-%d"),
              entry.status = entry_status_factor,
              exit.status = exit_status_factor,
              id=SID,
              data=cohort_df)
# NOTE: Dropping  249  rows with duration of follow up < tol

summary(lexis_df)
timeScales(lexis_df)

lexis_split_df<-cutLexis(data = lexis_df,
                  cut=cal.yr(lexis_df$ostodate,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK",
                  new.scale = "AK.Start"
)
summary(lexis_split_df)


lexis_df2<-cutLexis(data = lexis_split_df,
                  cut=cal.yr(lexis_split_df$LastAKdateplus120days,format = "%Y-%m-%d"),
                  timescale = "per",
                  new.state = "AK.quitted",
                  new.scale = "AK.quit"
)

summary(lexis_df2)

library(cohorttools)

boxesLx(lexis_df2,show.persons=FALSE)


# Split by calendar year (per)
range(lexis_df2$per)
lexis_df3<-splitLexis(lex=lexis_df2,time.scale = "per",breaks = c(2009,2011,2013,2015,2017))


# Year as factor
lexis_df3$per.c<-timeBand(lex = lexis_df3,time.scale = "per",type="factor")
# Year numeric
lexis_df3$per.num<-with(lexis_df3,per+lex.dur/2)

# Age as factor, age in middle of time slice
apu<-with(lexis_df3,cut(age+lex.dur/2,c(0,40,50,60,70,80,90,100,Inf)))
lexis_df3$age.c<-apu
table(apu,lexis_df3$lex.Xst)

# Muuttujien nimet
names(lexis_df3)

# Sex
apu<-factor(unclass(lexis_df3$SukupuoliBin),levels=0:1,labels =c("female","male"))
lexis_df3$sex<-apu

lexis_df3$sex <- relevel(lexis_df3$sex, ref = "male")

# Other variables
apu<-unclass(lexis_df3$HyperlipidemiaBOAC);table(apu)
lexis_df3$Hyperlipidemia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HypertensionBOAC);table(apu)
lexis_df3$Hypertension<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DiabetesBOAC);table(apu)
lexis_df3$Diabetes<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$HeartFailureBOAC);table(apu)
lexis_df3$HF<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalLiverFunctionBOAC);table(apu)
lexis_df3$LiverVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AbnormalRenalFunctionBOAC);table(apu)
lexis_df3$RenalVT<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$BleedingsBOAC);table(apu)
lexis_df3$Bleeding<-factor(apu,levels=0:1,labels =c("no","yes"))



apu<-unclass(lexis_df3$AlcoholBOAC);table(apu)
lexis_df3$ALKO<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$AnyVascularDiseaseBOAC);table(apu)
lexis_df3$anyvasc<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$CancerBeforeOrAtCohort);table(apu)
lexis_df3$syopa<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$DementiaBOAC);table(apu)
lexis_df3$dementia<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$PsychiatricDiseaseBOAC);table(apu)
lexis_df3$psykiatria<-factor(apu,levels=0:1,labels =c("no","yes"))

apu<-unclass(lexis_df3$tulot3luokkaa);table(apu)
lexis_df3$tulot3luokkaa<-factor(apu,levels=1:3,labels =c("low","mid", "high"))

apu<-unclass(lexis_df3$IschemicStrokeBOAC);table(apu)
lexis_df3$stroke<-factor(apu,levels=0:1,labels =c("no","yes"))



library(cohorttools)

#------------------------------
# Rate table
#------------------------------

tmp.rt<-mkratetable(Surv(lex.dur,lex.Xst=="AnyBleed")~per.c+age.c+Relevel(lex.Cst,list(1,2,3,4))+
                      sex+Hyperlipidemia+Hypertension+Diabetes+psykiatria+dementia+syopa+anyvasc+DM1onlyins2insoral3oral4noMed+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
                    data=lexis_df3,scale=100,add.RR = TRUE)

sink(file="RateTableUnique_Drugs_ATC_CategoryAnyBleedDRUGSAFTER.html")
knitr::kable(tmp.rt,caption="Rate table (1/100 person years)",format="html",digits=3)
sink()

#------------------------------
# Poisson models
#------------------------------

tmp.m1 <- glm(cbind(lex.Xst=="AnyBleed",lex.dur) ~ ALKO+Diabetes+per.c+Relevel(lex.Cst,list(1,2,3,4))+sex+age.c+tulot3luokkaa+Hypertension+Hyperlipidemia+psykiatria+dementia+syopa+anyvasc+Unique_Drugs_ATC_Category+Bleeding+stroke+LiverVT+RenalVT+HF,
              data=lexis_df3,family="poisreg")


round(ci.exp(tmp.m1),2)

tmp.m2 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:sex)
tmp.m3 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:tulot3luokkaa)
tmp.m4 <- update(tmp.m1, ~ .+Unique_Drugs_ATC_Category:Relevel(lex.Cst,list(1,2,3,4)))

SEXinteractiomAnyBleed <- anova(tmp.m2,tmp.m1,test="Chisq")
TULOinteractiomAnyBleed <- anova(tmp.m3,tmp.m1,test="Chisq")
AKinteractiomAnyBleed <- anova(tmp.m4,tmp.m1,test="Chisq")

round(ci.exp(tmp.m2),2)
round(ci.exp(tmp.m3),2)
round(ci.exp(tmp.m4),2)

sink(file="IRRTulosUnique_Drugs_ATC_CategoryAnyBleedDRUGSAFTER.html")

knitr::kable(round(ci.exp(tmp.m1),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m2),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m3),2),caption="IRR from Poisson regression model",format="html")
knitr::kable(round(ci.exp(tmp.m4),2),caption="IRR from Poisson regression model",format="html")

sink()

