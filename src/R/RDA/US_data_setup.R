library(stringr)
library(magrittr)
library(spdep)
library(maps)
library(maptools)
library(data.table)

age_adjusted <- fread(file.path(getwd(), "data", "US_data", "age_adjusted.csv"))
colnames(age_adjusted) <- c("state_county", "site", "age_adjusted_rate", "count", "population")
age_adjusted[, uniqueN(state_county)]
lung <- age_adjusted[site == "Lung and Bronchus",
                     .(state = str_match(state_county, "^(\\w{2}).+")[,2],
                       county_name = str_match(state_county, ":\\s(.+)\\s\\(\\d{5}\\)")[,2],
                       county_fips = str_match(state_county, ":\\s.+\\s\\((\\d{5})\\)")[,2] %>%
                         as.integer(),
                       age_adjusted_rate, count, population)]
lung <- lung[!is.na(county_name)]
# Import US counties
county_poly <- maps::map("county", fill = TRUE, plot = FALSE)
county_state <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[1]]))
county_names <- strsplit(county_poly$names, ",") %>%
  sapply(function(x) str_to_title(x[[2]]))
sf_use_s2(FALSE)
county_df <- data.table(state_name = county_state,
                        state_abb = state.abb[match(county_state, state.name)],
                        county_fips = county.fips[match(county_poly$names, county.fips$polyname),]$fips,
                        county_name = county_names)
county_sp <- maptools::map2SpatialPolygons(county_poly, IDs = county_poly$names)
county_nbs <- poly2nb(county_sp)
no_neighbors <- vapply(county_nbs, function(x) identical(x, 0L), logical(1))
# restrict to connected county map
county_sp <- county_sp[!no_neighbors,]
county_state <- county_state[!no_neighbors]
county_names <- county_names[!no_neighbors]
county_nbs <- poly2nb(county_sp)
W <- nb2mat(county_nbs, style = "B")
D <- diag(rowSums(W))
alpha <- 0.99
Q <- D
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
scaling_factor <- exp(mean(log(diag(Sigma))))
Q_scaled <- Q * scaling_factor
Q_scaled_cholR <- Q_cholR * sqrt(scaling_factor)
Sigma_scaled <- Sigma / scaling_factor
N <- nrow(W)

# data only counties subset?

data <- merge(county_df[!no_neighbors], lung, by = c("county_fips"), all.x = TRUE)
has_data <- cbind(has_data = data[, !is.na(population)], st_as_sf(county_sp))
plot(has_data,
     main = paste0("Data In: ", paste0(lung[, unique(state)], collapse = ", ")))
