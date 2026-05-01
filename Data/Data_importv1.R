library(dplyr)
library(ggplot2)

# fishery (large-purse-seine) catch
#	Yield (MT), CPUE (MT/haul)
data_yield_mackerel_all = read.csv("D:\\Project\\Thysis\\Data\\Age_structured_data\\Yearly_yield_by_region.csv", header=T);  # Chub mackerel yield by region
data_yield_mackerel = as.matrix(data_yield_mackerel_all);
data_cpue_mackerel_all = read.csv("D:\\Project\\Thysis\\Data\\Age_structured_data\\Yearly_CPUE_by_region.csv", header=T);  # Chub mackerel CPUE by region
data_cpue_mackerel <- as.matrix(data_cpue_mackerel_all);

Total_yield=data_yield_mackerel[,"Total"];             # vector, Yt: Yield in MT at year t
East_yield=data_yield_mackerel[,"East"];
West_yield=data_yield_mackerel[,"West"];
South_yield=data_yield_mackerel[,"South"];

Total_CPUE=data_cpue_mackerel[,"Total"]
East_CPUE=data_cpue_mackerel[,"East"]
West_CPUE=data_cpue_mackerel[,"West"]
South_CPUE=data_cpue_mackerel[,"South"]

# number of age classes
nages=6;   

# number of years for the fishery catch data, cpue data, and length frequency data: nyrs of 2000-2024;
nyrs=dim(data_yield_mackerel)[1]

# length frequency data (23 x 43);
# 25: year 2000 - 2024
# 43: midpoints (cm) of length classes: 10.5, 11.5, 12.5, ... , 52.5;
data_length_freq_annual<-read.csv("D:\\Project\\Data\\lengthbased_data\\data_length_frequency_data_250407.csv", header=T);
data_total_freq_annual<-rowSums(data_length_freq_annual)
# data_length_freq_annual<-as.matrix(data_length_freq_annual); #including data from 2000
# data_length_freq_annual; ### matrix
# data_length_prob_mackerel<-data_length_freq_mackerel/apply(data_length_freq_mackerel, 1, sum)
# data_length_freq_first_half<-read.csv("D:\\Project\\Data\\lengthbased_data\\data_length_freq_first_half_250407.csv", header=T);
# data_length_freq_first_half<-as.matrix(data_length_freq_first_half); #include data from 2000
# data_length_freq_first_half
# data_length_freq_second_half<-read.csv("D:\\Project\\Data\\lengthbased_data\\data_length_freq_second_half_250407.csv", header=T);
# data_length_freq_second_half<-as.matrix(data_length_freq_second_half); #include data from 2000
# data_length_freq_second_half

# Age composition data (25 x 6)
data_age_comp_annual_all=read.csv("D:\\Project\\Thysis\\Data\\Age_structured_data\\Catch_at_age_using_FIALkey.csv", header=T);
data_age_comp_annual=data_age_comp_annual_all[,-2] #Age 1-6+
data_age_comp_annual

## Prior info
# Korean mackerel
#bounds_q: lower and upper bounds of q for "logit( (q-lower)/(upper-lower))"
bounds_q=c(exp(-40.0), exp(-10.0));

# cm
mean_length_chub <- c(27.8, 30.3, 34.8, 37.2, 40.1, 42.7)

# Weight
alpha_LW_Total <- 0.00269; beta_LW_Total <- 3.453
alpha_LW_West  <- 0.00273; beta_LW_West  <- 3.455
alpha_LW_South <- 0.00245; beta_LW_South <- 3.479
alpha_LW_East  <- 0.00427; beta_LW_East  <- 3.322

# Mean weight at age a in area y
mean_WAA_Total <- (alpha_LW_Total * mean_length_chub ^ beta_LW_Total) / 1000
mean_WAA_West  <- (alpha_LW_West  * mean_length_chub ^ beta_LW_West)  / 1000
mean_WAA_South <- (alpha_LW_South * mean_length_chub ^ beta_LW_South) / 1000
mean_WAA_East  <- (alpha_LW_East  * mean_length_chub ^ beta_LW_East)  / 1000

# Print results
print(round(mean_WAA_Total, 4))
print(round(mean_WAA_West,  4))
print(round(mean_WAA_South, 4))
print(round(mean_WAA_East,  4))

# data_LW=read.csv(file = "D:\\Project\\Thysis\\Data\\Age_structured_data\\length_weight_data.csv", header = T);
# data_LW<-as.matrix(data_LW);
# data_LW; # matrix
# 
# data_length_LW=data_LW[,"length.cm."]                   ### vector, length_data: Fork length data (cm)
# data_weight_LW=data_LW[,"weight.g."]                    ### vector, weigth_data: weight data (gram)

# plot(data_length_LW, data_weight_LW, xlab="Fork length (cm)", ylab="Weigth (g)", main="Length-weight data", cex.lab=1.5, cex.axis=1.5, cex.main=2.0)

# log_alpha_LW=log(0.0028);
# log_beta_LW=log(3.43);
# 
# mean_length_chub=c(27.8,30.3,34.8,37.2,40.1,42.7)
# mean_WAA=0.0028*mean_length_chub^3.43/1000
# mean_WAA   #gram

# length(cm)-maturation(ratio) data (scanned from the paper of Sora Kim et al. 2020)
#data_maturation=read.csv(file = "maturation_data.csv", header = T);
#data_maturation<-as.matrix(data_maturation);
#data_maturation; ### matrix
# data=read.csv(file = "D:\\Project\\Thysis\\Data\\Age_structured_data\\Maturation_data_20250420.csv", header = T);  # 1 - 12
# 
# names(data)=c("Year","Month","FL","BW","Mat")
# 
# data <- data %>%
#   mutate(Mat_binary = ifelse(Mat %in% 1, 0, 1))    # case 1
# data <- data %>% filter(Month >= 1 & Month <= 6)
# 
# breaks <- seq(10, 53, by = 1)
# data$FL_class <- cut(data$FL, breaks = breaks, right = FALSE)
# 
# midpoints <- head(breaks, -1) + 0.5  
# levels(data$FL_class) <- midpoints  
# data$FL_mid <- as.numeric(as.character(data$FL_class))  
# 
# logistic_data <- data %>%
#   group_by(FL_mid) %>%
#   summarise(
#     m = n(),           
#     y = sum(Mat_binary)      
#   )
# 
# full_logistic_data <- data.frame(FL_mid = bin) %>%
#   left_join(logistic_data, by = "FL_mid") %>%
#   mutate(
#     m = ifelse(is.na(m), 0, m),
#     y = ifelse(is.na(y), 0, y)
#   )
# 
# data_maturation=as.matrix(full_logistic_data)
# data_maturation
# 
# #data_length_maturation=data_maturation[,"Fork_length.cm."]    ### vector, mat_length_data: maturation fork length (cm)
# #data_rate_maturation=data_maturation[,"maturation_rate"]      ### vector, mat_ratio_data: maturation rate (0 ~ 1)
# 
# # plot(data_length_maturation, data_rate_maturation, xlab="Fork length (cm)", ylab="Rate (0.0 ~ 1.0)", main="Length-maturation data", cex.lab=1.5, cex.axis=1.5, cex.main=2.0)
# 
# b0_mat=20.11;
# b1_mat=0.7;
# #b0_mat=8.695706
# #b1_mat=0.244470 
# x=seq(10, 53, 1)
# mature=1/(1+exp(b0_mat-b1_mat*x))
# #plot(x,mature)
# # ratio_female ####
# ratio_female=0.6;

# length(cm)-fecundity(the number of eggs) data (achieved by sanning data in Cha et al. 2002) ####
data_fecundity=read.csv(file = "D:\\Project\\Thysis\\Data\\Age_structured_data\\fecundity_data.csv", header = T);
data_fecundity<-as.matrix(data_fecundity);
data_fecundity; ### matrix

data_length_fecundity=data_fecundity[,"Fork_length.cm."]     ### vector, fec_length_data: fecundity fork length data (cm)
data_eggs_fecundity=data_fecundity[,"the.number.of.eggs"]    ### vector, fec_egg_data: fecundity rate data (the number of eggs)

log_fec_a=-0.015;
log_fec_b=1.26;


####################################################################
####################################################################
####################################################################
# dome-shape selectivity by age ===> dome-shape selectivity by length

load("D:\\Project\\Thysis\\Data\\Age_structured_data\\kormackerel.RData");

#JW's study;
select_dome_estimates =c(1.392, 5.349, 0.355, 0.661);  #estimates of four parameters in the dome-selectivity;

inf1=select_dome_estimates[1];
inf2=select_dome_estimates[2];
g1=select_dome_estimates[3];
g2=select_dome_estimates[4];
sel_dome=c();
for(age in 1:6)
  sel_dome[age]=(1/(1+exp(-age+inf1)/g1) )*(1/(1+exp(age-inf2)/g2) );
