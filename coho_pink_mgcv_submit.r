# First version of GAM modeling for pink salmon effects on coho
# This version is Salish Sea only
library("tidyverse")
library("mgcv")
library("gratia")

setwd("C:/Users/nelso/OneDrive/BNelson/Consulting/Neala WDFW/2020 project")

## Clear everything ##
rm(list=ls()) 
gc() 
# Set randoms #
set.seed(99)

##-------------------------------------------------------------------------------------------------------------------------------------------------------------
# Global functions
##-------------------------------------------------------------------------------------------------------------------------------------------------------------
source("pink_analysis_functions.r")
##-------------------------------------------------------------------------------------------------------------------------------------------------------------
# Import Data
##-------------------------------------------------------------------------------------------------------------------------------------------------------------
survival_data<- read.csv("coho_survival_updated_feb_2021.csv", header=TRUE, stringsAsFactors = FALSE)
release_data<- read.csv("coho_releases_updated.csv", header = TRUE, stringsAsFactors = FALSE)
pink_abundance<- read.csv("pink_abundance_2021.csv", header=TRUE)
pink_fry<- read.csv("salish_pink_fry.csv", header=TRUE)
covariates<- read.csv("additional_covariates.csv", header=TRUE)
zooplankton<- read.csv("zooplankton_indices.csv", header=TRUE)

## Covariates pre-processing
covariates$salish_seals<- log(covariates$sog_seals + covariates$ps_seals + covariates$juan_seals) # Salish Sea seals
covariates$coast_seals<- log(covariates$coast_seals) # Coastal seals
covariates$columbia_discharge<- (covariates$columbia_discharge/100000) # Columbia R discharge
covariates$salish_herring<- log(covariates$salish_herring) # Salish Sea herring
covariates$sst<- covariates$sst - mean(covariates$sst, na.rm = TRUE)
pink_fry$ps<- pink_fry$ps/1000000 # Convert to millions of fry
pink_fry$year<- pink_fry$year + 1 # Change to migration year
pink_fry$total_fry<- pink_fry$ps + pink_fry$fraser # Add col for total Salish Sea fry (Fraser + PS)
pink_fry<- na.omit(pink_fry)

##-------------------------------------------------------------------------------------------------------------------------------------------------------------
# Lists and labels
##-------------------------------------------------------------------------------------------------------------------------------------------------------------
all_regions<- c("JUAN", "HOOD", "SPS", "MPS", "NPS", "NOWA", "SOG", "FRA", "COL", "WILP", "WCVI", "NWC")
ps_regions<- c("JUAN", "HOOD", "SPS", "MPS", "NPS", "NOWA")
sog_regions<- c("SOG", "FRA")
salish_regions<- c(ps_regions, "SOG", "FRA")
coastal_regions<- c("WILP", "NWC", "WCVI")
columbia_regions<- c("COL")

# Remove excluded stocks
survival_data<- survival_data %>% filter(exclude=="No")
## Modify survival data
survival_data<- na.omit(survival_data) # Get rid of entries (rows) with NAs
# Add 'release year' col
survival_data$release_year<- survival_data$year
# Set first and last years in the release dataset that have both years so the following loop will execute successfully (i.e. no NAs)
min_year<- min(survival_data$release_year)
max_year<- max(survival_data$release_year)
# Constrain years in analysis
survival_data<- filter(survival_data, release_year>=min_year, release_year<=max_year)
# List of stock IDs
stock_id<- unique(survival_data$stock_id)
full_stock_names<- unique(survival_data$stock_name)
# Even years are when pink salmon are present
survival_data$pink_year<- ifelse(survival_data$release_year%%2==0, 1, 0)
# Add region code
survival_data$region_code<- NA
regions<- unique(survival_data$region_id)
for(i in 1:length(regions)) {survival_data[survival_data$region_id==regions[i], "region_code"]<- i } 
# Add population code
survival_data$stock_id<- as.character(survival_data$stock_id)  
ps_pops<- unique(survival_data$stock_id)
survival_data$stock_code<- NA
for(i in 1:length(ps_pops)) {survival_data[survival_data$stock_id==ps_pops[i], "stock_code"]<- i } 
# Add year effect 
survival_data$year_effect<- survival_data$release_year-min(survival_data$release_year)+1
# Create the release values
release_data<- filter(release_data, release_year >= min_year & release_year <= max_year)
release_data<- release_data %>% 
  group_by(region) %>%
  mutate(z_releases=sd_std(releases, releases))
survival_data<- left_join(survival_data, release_data, by=c("release_year", "basin"="region")) %>%
  select(-releases)
# Add abundance of pinks
pink_year_seq<- data.frame(year=seq(from=(min(pink_fry$year)-1), to=(max(pink_fry$year)+1), by=1))
pink_year_seq<- full_join(pink_year_seq, pink_fry, by="year")
pink_year_seq[is.na(pink_year_seq$total_fry), "total_fry"]<- 1
survival_data<- left_join(survival_data, pink_year_seq, by=c("release_year"="year"))
# Calculate instantaneous mort.
survival_data$instant_m<- -1*log(survival_data$surv)
# Select dataset for analysis
survival_data<- filter(survival_data, dataset_id==1)
#survival_data<- filter(survival_data, region_id %in% columbia_regions)
# Number of observations/data points
n_obs<- nrow(survival_data)

##-------------------------------------------------------------------------------------------------------------------------------------------------------------
# Build data matrix for mgcv model-building
##-------------------------------------------------------------------------------------------------------------------------------------------------------------
model_data<- data.frame(instant_m=as.numeric(survival_data$instant_m),
                        pinks=as.factor(survival_data$pink_year),
                        pink_abundance=log(as.numeric(survival_data$total_fry)),
                        releases=as.numeric(survival_data$z_releases),
                        basin=as.factor(survival_data$basin),
                        region=as.factor(survival_data$region_code),
                        stock=as.factor(survival_data$stock_id),
                        stock_name=as.character(survival_data$stock_name),
                        origin=as.factor(survival_data$origin),
                        year_code=as.factor(survival_data$year_effect),
                        year=as.numeric(survival_data$release_year),
                        year_dig=substr(as.character(as.numeric(survival_data$release_year)), 3, 4),
                        stock_lab=survival_data$stock_name)

# Add additional covariates
covariates<- covariates %>% select(-quazi_bi_std, -columbia_discharge)
model_data<- left_join(model_data, covariates, by="year") %>%
  drop_na() %>%
  mutate(releases_x_pinks=(as.numeric(pinks)-1)*releases, releases_x_pink_ab=pink_abundance*releases)

## Summarize stock data
stock_data<- NULL
for(i in 1:length(unique(model_data$stock))){
  temp_stock<- filter(model_data, stock_name==full_stock_names[i])
  temp_n<- nrow(temp_stock)
  temp_year_range<- paste(min(temp_stock$year), "-", max(temp_stock$year), " (", temp_n, ")", sep="")
  temp_origin<- as.character(unique(temp_stock$origin))
  temp_basin<- as.character(unique(temp_stock$basin))
  temp_name<- full_stock_names[i]
  temp_outs<- c(temp_name, temp_basin, temp_origin, temp_year_range)
  if(i==1){stock_data<- temp_outs} else{
    stock_data<- rbind(stock_data, temp_outs)
  }
}

# Build stock df
stock_data<- as.data.frame(stock_data); names(stock_data)<- c("Stock name", "Region", "Origin", "Years")
stock_data<- arrange(stock_data, desc(Region))
#write.csv(stock_data, "stock_data_coho.csv", row.names = FALSE)

## Model fitting section
# Salish Sea
salish_data<- model_data[model_data$basin %in% c("PS", "SOG"), ]
# Linear
mod_0<- gamm(log(instant_m) ~ origin + salish_seals + sst + salish_herring + npgo + pinks + releases + releases_x_pinks + s(stock, k=length(unique(salish_data$stock)), bs="re"),
             data = salish_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_0$gam)
AIC(mod_0$lme)
draw(mod_0$gam)

# Linear
mod_0_1<- gamm(log(instant_m) ~ origin + salish_seals + sst + salish_herring + npgo + pink_abundance + releases + releases_x_pink_ab + s(stock, k=length(unique(salish_data$stock)), bs="re"),
               data = salish_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_0_1$gam)
AIC(mod_0_1$lme)
draw(mod_0_1$gam)

# GAM
mod_1<- gamm(log(instant_m) ~ origin + salish_seals + sst + salish_herring + npgo + s(releases, by=pinks, bs="tp") + s(stock, k=length(unique(salish_data$stock)), bs="re"),
             data = salish_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_1$gam)
AIC(mod_1$lme)
draw(mod_1$gam)

# GAM w/ interaction (pink years)
mod_1_2<- gamm(log(instant_m) ~ origin + salish_seals + sst + salish_herring + npgo + pinks + releases + te(year, releases_x_pinks, bs="tp") + s(stock, k=length(unique(salish_data$stock)), bs="re"),
               data = salish_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_1_2$gam)
AIC(mod_1_2$lme)
draw(mod_1_2$gam, select = 1, scales = "fixed")
int_term<- smooth_estimates(mod_1_2$gam) %>%
  group_by(year) %>%
  summarise(mean_est=mean(est), up_est=mean(est+(2*se)), low_est=mean(est-(2*se)))
# Plot interaction vs time

plot(int_term$year, int_term$mean_est, ann = FALSE, las=1, type="l", lwd=2, ylim=c(-2,1), col="purple")
abline(h=0.0, lty="dashed")
lines(int_term$year, int_term$low_est, lty="dashed", lwd=2, col="purple")
lines(int_term$year, int_term$up_est, lty="dashed", lwd=2, col="purple")
mtext(side=1, line=2.5, font=2, "Release year", cex=1.25)
mtext(side=2, line=2.75, font=2, "Pinks (even year) x releases", cex=1.25)

# GAM w/ interaction (pink years)
mod_1_3<- gamm(log(instant_m) ~ origin + salish_seals + sst + salish_herring + npgo + pink_abundance + releases + te(year, releases_x_pink_ab, bs="tp") + s(stock, k=length(unique(salish_data$stock)), bs="re"),
               data = salish_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_1_3$gam)
AIC(mod_1_3$lme)
draw(mod_1_3$gam, select = 1, scales = "fixed")

# Coastal
coastal_data<- model_data[model_data$basin=="COAST", ]

mod_2<- gamm(log(instant_m) ~ coast_seals + sst_arc + npgo + pinks + releases + releases_x_pinks + s(stock, k=length(unique(coastal_data$stock)), bs="re"),
             data = coastal_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_2$gam)
AIC(mod_2$lme)
draw(mod_2$gam)

mod_2_1<- gamm(log(instant_m) ~ coast_seals + sst_arc + npgo + pink_abundance + releases + releases_x_pink_ab + s(stock, k=length(unique(coastal_data$stock)), bs="re"),
               data = coastal_data, method = "REML", family=gaussian(link="identity"), correlation = corAR1(form = ~ year|stock))
summary(mod_2_1$gam)
AIC(mod_2_1$lme)
draw(mod_2_1$gam)

# Plot stocks
mod_1<- mod_1_2
preds<- predict(mod_1$gam, se.fit=TRUE)
salish_data <- transform(salish_data, preds_med = preds$fit, preds_ci = preds$se.fit)
theme_set(theme_bw()) 
theme_update(panel.grid = element_blank())
ggplot(data=salish_data, aes(x=year, y=instant_m, group=stock)) +
  facet_wrap(~stock_lab) +
  geom_ribbon(aes(ymin=exp(preds_med-2*preds_ci),
                  ymax=exp(preds_med+2*preds_ci)), alpha=0.25) +
  geom_line(aes(y=exp(preds_med))) +
  geom_point() +
  #facet_rep_grid(~stock, labeller = labeller(stock=stock_lab)) +
  labs(x="Release year", y="Instant. Mortality")

# Obs vs. pred
ggplot(salish_data, aes(x=exp(preds_med), y=instant_m)) +
  facet_wrap(~stock_lab) +
  geom_point(alpha=0.1) +
  geom_abline() +
  labs(x="Predicted M", y="Observed M")
