library(stringr)
library(magrittr)
library(spdep)
library(maps)
library(maptools)
library(data.table)
library(openxlsx)

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

# Fix D.C.
county_df[county_fips == 11001, `:=`(state_name = "Washington",
                                     county_name = "District Of Columbia")]

# Shannon County, SD (FIPS code = 46113) was renamed Oglala Lakota County and
# assigned a new FIPS code (46102) effective in 2014.
county_df[county_name == "Oglala Lakota", county_fips := 46102]
county_sp <- maptools::map2SpatialPolygons(county_poly, IDs = county_df$county_fips)
setkey(county_df, county_fips)
county_sf <- st_as_sf(county_sp)
st_crs(county_sf) <- st_crs(st_as_sf(county_poly))
county_sf$fips <- unique(county_df$county_fips)
county_nbs <- poly2nb(county_sp)
no_neighbors <- vapply(county_nbs, function(x) identical(x, 0L), logical(1))

# restrict to connected county map
county_sp <- county_sp[!no_neighbors,]
county_state <- county_state[!no_neighbors]
county_names <- county_names[!no_neighbors]
county_sf <- county_sf[!no_neighbors,]

# cancer and smoking data
cancer <- read.xlsx(file.path(getwd(), "data", "US_data", "IHME_county_cancer_mortality.xlsx"),
                    sheet = "Tracheal, bronchus, & lung ", startRow = 2)
mortality_col_names <- paste0("mortality", c(1980, 1985, 1990, 1995, 2000, 2005,
                                             2010, 2014))
colnames(cancer) <- c("location", "county_fips", mortality_col_names, "total_change")
setDT(cancer)
# take out confidence intervals
confintCols <- c(mortality_col_names, "total_change")
extractRate <- function(x) as.numeric(str_match(x, "^(.+)\\s\\(")[,2])
cancer[, (confintCols) := lapply(.SD, extractRate),
       .SDcols = c(confintCols)]
cancer[, county_fips := as.integer(county_fips)]
cancer[, is_county := county_fips >= 1000]
cancer <- cancer[!is.na(county_fips) & is_county & county_fips %in% county_sf$fips,]

has_data_sf <- cbind(county_sf, has_data = county_sf$fips %in% cancer$county_fips)
#plot(has_data_sf[, "has_data"])

smoking <- fread(file.path(getwd(), "data", "US_data", "IHME_data",
                           "IHME_US_COUNTY_TOTAL_AND_DAILY_SMOKING_PREVALENCE_1996_2012.csv"))
smoking <- smoking[year == 2012 & sex == "Both" & county != "",
                   .(state, county, total_mean_smoking = total_mean,
                     daily_mean_smoking = daily_mean,
                     location = paste(county, state, sep = ", "))]
smoking[county == "Shannon County" & state == "South Dakota",
        location := "Oglala County, South Dakota"]
smoking[, location := str_replace(location, "St.", "Saint")]
# cancer[!location %in% smoking$location]
cancer_smoking <- merge(cancer, smoking, by = "location")
cancer_smoking$is_county <- NULL
cancer_smoking_sf <- merge(county_sf, cancer_smoking,
                           by.x = "fips", by.y = "county_fips")
cancer_smoking <- st_drop_geometry(cancer_smoking_sf)
smoking[, county_name := str_match(county, "^(.+)\\s(County|City|Parish)")[,2]]
has_data_sf <- cbind(has_data = county_sf$fips %in% cancer_smoking_sf$fips, county_sf)
# plot(has_data_sf[, "has_data"], main = "Lung Cancer Mortality and Smoking Data Coverage")

# absolutely make sure data regions are connected by subsetting to largest piece
# tol <- .01
# buffered <- st_buffer(cancer_smoking_sf, tol)
# all_geoms <- st_union(buffered)
# all_geoms <- st_cast(all_geoms, "POLYGON")
# largest_part <- all_geoms[which.max(st_area(all_geoms))]
# contig_indx <- st_intersects(cancer_smoking_sf, largest_part, sparse = FALSE)
# contig_data <- cancer_smoking_sf[contig_indx, ]
# plot(buffered[contig_indx, "fips"])

county_nbs <- poly2nb(cancer_smoking_sf)
W <- nb2mat(county_nbs, style = "B")
D <- diag(rowSums(W))
alpha <- 0.99
Q <- D
Q <- D - alpha * W
Q_cholR <- chol(Q)
Sigma <- chol2inv(Q_cholR)
scaling_factor <- exp(mean(log(diag(Sigma))))
Q_scaled <- Q * scaling_factor
N <- nrow(W)

ij_list <- data.frame(
  i = rep(seq_len(N), times = vapply(county_nbs, length, numeric(1))),
  j = unlist(county_nbs)
)
ij_list <- ij_list[ij_list$i < ij_list$j, ]
rownames(ij_list) <- NULL

rm(cancer, smoking, county_poly,
   county_names, county_state, county_sp, county_df,
   confintCols, mortality_col_names, no_neighbors, extractRate)
rm(W, D, alpha, Sigma, Q_cholR, Q)
gc()
