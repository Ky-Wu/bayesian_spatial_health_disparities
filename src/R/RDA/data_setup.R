library(data.table)
library(stringr)
library(magrittr)
library(spdep)
library(maps)
library(maptools)

covariates <- fread(file.path(getwd(), "data", "covariates.csv"))
race <- fread(file.path(getwd(), "data", "race.csv"))
sex <- fread(file.path(getwd(), "data", "sex.csv"))
insurance <- fread(file.path(getwd(), "data", "insurance.csv"))
smoking <- fread(file.path(getwd(), "data", "smoking.csv"))
SIR_adjusted <- fread(file.path(getwd(), "data", "SIR_adjusted.csv"))

covariates <- covariates[, .(
  state = str_match(State_county, "^(\\w{2}).+")[,2],
  county_name = str_match(State_county, ":\\s(.+)\\s(County|Registry)\\s\\(")[,2],
  county_fips = str_match(State_county, ":\\s.+\\s\\((\\d{5})\\)")[,2],
  young = V_Persons_age_18_ACS_2012_2016 / 100,
  old = V_Persons_age_65_ACS_2012_2016 / 100,
  highschool = VHighschooleducationACS2012201 / 100,
  poverty = VFamiliesbelowpovertyACS201220 / 100,
  unemployed = V_Unemployed_ACS_2012_2016 / 100
)]
covariates <- covariates[state == "CA"]
black_percents <- race[
  substr(race$State_county, 1, 2) == "CA" &
    race$Race_recode_White_Black_Other=="Black",
  .(state = str_match(State_county, "^(\\w{2}).+")[,2],
    county_name = str_match(State_county, ":\\s(.+)\\s(County|Registry)\\s\\(")[,2],
    county_fips = str_match(State_county, ":\\s.+\\s\\((\\d{5})\\)")[,2],
    black_percent = Row_Percent)]
male_percents <- sex[
  substr(sex$State_county, 1, 2) == "CA" &
    sex$Sex == "Male",
  .(state = str_match(State_county, "^(\\w{2}).+")[,2],
    county_name = str_match(State_county, ":\\s(.+)\\s(County|Registry)\\s\\(")[,2],
    county_fips = str_match(State_county, ":\\s.+\\s\\((\\d{5})\\)")[,2],
    male_percent = Row_Percent)]
uninsured_percents <- insurance[
  substr(insurance$State_county, 1, 2) == "CA" &
    insurance$Insurance_Recode_2007 == "Uninsured",
  .(state = str_match(State_county, "^(\\w{2}).+")[,2],
    county_name = str_match(State_county, ":\\s(.+)\\s(County|Registry)\\s\\(")[,2],
    county_fips = str_match(State_county, ":\\s.+\\s\\((\\d{5})\\)")[,2],
    uninsured_percent = Row_Percent)]
colnames(smoking) <- c("county_name", "smoking")
smoking[, `:=`(smoking = str_match(smoking, "(.+)%")[,2],
               county_name = str_to_title(str_trim(county_name)))]
covariates_all <- covariates %>%
  merge(black_percents, by = c("state", "county_name", "county_fips")) %>%
  merge(male_percents, by = c("state", "county_name", "county_fips")) %>%
  merge(uninsured_percents, by = c("state", "county_name", "county_fips")) %>%
  merge(smoking, by = c("county_name"))

covariate_columns <- c("young", "old", "highschool", "poverty", "unemployed",
                       "black_percent", "male_percent", "uninsured_percent",
                       "smoking")
covariates_all[, (covariate_columns) := lapply(.SD, as.numeric),
               .SDcols = covariate_columns]
lung <- SIR_adjusted[Site.code == "Lung and Bronchus",
                     .(state = str_match(State.county, "^(\\w{2}).+")[,2],
                       county_name = str_match(State.county, ":\\s(.+)\\s(County|Registry)\\s\\(")[,2],
                       county_fips = str_match(State.county, ":\\s.+\\s\\((\\d{5})\\)")[,2],
                       O_count, E_count, standard_ratio)]
lung <- merge(lung, covariates_all, by = c("state", "county_name", "county_fips"))

# Read in CA county map and compute adjacency matrix W and CAR prec. matrix Q
county_poly <- maps::map("county","california", fill=TRUE, plot=FALSE)
county_names <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[2]]))
sf_use_s2(FALSE)
county_sf <- st_as_sf(county_poly)
rownames(county_sf) <- NULL
county_sf$county_name <- county_names
county_sf <- merge(county_sf, lung, by = c("county_name"))
county_nbs <- poly2nb(county_sf, queen = FALSE)
W <- nb2mat(county_nbs, style="B")
D <- diag(rowSums(W))
alpha <- 0.99
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
scaling_factor <- exp(mean(log(diag(Sigma))))
Q_scaled <- Q * scaling_factor
Q_scaled_cholR <- Q_cholR * sqrt(scaling_factor)
Sigma_scaled <- Sigma / scaling_factor
N <- nrow(W)

nbs_df <- data.frame(
  node1 = rep(seq_len(N), times = lapply(county_nbs, length)),
  node2 = unlist(county_nbs)
)
ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL
