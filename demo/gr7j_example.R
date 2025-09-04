# Trying to get GR4J rainfall runoff model working


# load package - not necessary to do each time when working in a project
library(hydromad)

# Load the new GR4J_rory model. This loads the altered functions into the environment
# source("GR4J_rory.R")
# hydromad.model("gr4j_rory", sma = gr4j_rory.sim)

hydromad.options(pure.R.code = TRUE)  # ensures we use rcode for gr4j, not C

# Load the HydroTestData dataset. Contains P, E and Q values for a timeseries 2000 Jan to March
data(HydroTestData)



# BUILD the model using the hydromad() function; calling the GR4J sma and routing functions
mod0 <- hydromad(
  HydroTestData, 
  sma = "gr4j_roryVar", #gr4j_roryVar
  routing = "gr4jrouting",
  x1 = 100,
  x2 = 20, x3 = 1, x4 = 10,
  x5 = 0.15, x6 = 0.4, x7 = 0.25
)

testQ <- predict(mod0, return_state = TRUE)

## plot results with state variables
xyplot(
  cbind(HydroTestData[, 1:2], gr4j = testQ),
  strip.left = TRUE,  # puts strip labels on the left side
  strip = strip.custom(
    factor.levels = c(
      "Precipitation (mm/day)",
      "Potential ET (mm/day)",
      "Effective Rainfall (U, mm/day)",
      "Soil Moisture (S, mm)",
      "Actual ET (mm/day)",
      "Routing Store (X, mm)"
    )
  ),
  scales = list(y = list(relation = "free")),
  xlab = "Date"
)


## ---- Plot AET/PET vs S/SMSC ----

# Extract time series and PET
time_vals <- index(testQ)
PET <- coredata(HydroTestData[, "E"])

# Find valid time steps: PET > 0 and no NA in ET or S
valid_idx <- which(
  !is.na(testQ$S) &
    !is.na(testQ$ET) &
    PET > 0
)

# Create a data frame for plotting
df_plot <- data.frame(
  S_over_SMSC = coredata(testQ$S[valid_idx]) / coef(mod0)["x1"],
  AET_over_PET = coredata(testQ$ET[valid_idx]) / PET[valid_idx],
  Time = as.Date(time_vals[valid_idx])
)

# Load ggplot2
library(ggplot2)

# ~~~~~~~START Note: Make curve to see if it matches produced AET/PET vs S/SMSC
# Generated vs observational
curve_df <- data.frame(S_over_SMSC = seq(0, 1, length.out = 500))

# Extract parameter values from the model
x5 <- coef(mod0)["x5"]
x6 <- coef(mod0)["x6"]
x7 <- coef(mod0)["x7"]

# Compute curve values
curve_df$AET_over_PET <- (curve_df$S_over_SMSC)^x5 - x6 * (1 - (1 - curve_df$S_over_SMSC)^x7)
# ~~~~~~~~END

# Plot using ggplot2 with date-based colour scale and overlay model curve
ggplot(df_plot, aes(x = S_over_SMSC, y = AET_over_PET, colour = Time)) +
  geom_point(size = 1.5) +
  geom_line(data = curve_df, aes(x = S_over_SMSC, y = AET_over_PET),
            inherit.aes = FALSE,
            colour = "black", linetype = 'dashed', alpha = 0.5, linewidth = 1.1) +
  scale_colour_viridis_c(
    option = "D",
    name = "Date",
    guide = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 1,
      barheight = 10
    ),
    labels = function(x) format(as.Date(x, origin = "1970-01-01"), "%d %b %Y"),
    breaks = pretty(as.numeric(df_plot$Time), n = 5)
  ) +
  labs(
    title = "AET/PET vs S/SMSC (Coloured by Date + ET Curve)",
    x = "S / SMSC",
    y = "AET / PET"
  ) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_minimal()

#######Pareto analyis and caliberation for the new model


library(hydromad)
library(readxl)

# Load SAFT drought function
source("C:/Users/thris/Desktop/FYP-B/saft_drought_algorithm.R", chdir=TRUE)

# Load data
catchment_data <- read_excel("C:/Users/thris/Desktop/FYP-B/Camels Data Extract.xlsx", sheet = 1)
rainfall_data <- as.vector(catchment_data$P)
annual_rainfall <- aggregate(P ~ year, data = catchment_data, sum)
is_drought <- saft_drought_algorithm(annual_rainfall$P)
dry_years <- annual_rainfall$year[is_drought]
nondry_years <- annual_rainfall$year[!is_drought]

# Load full P-E-Q timeseries
catchment_data_sheet2 <- read_excel("C:/Users/thris/Desktop/FYP-B/Camels Data Extract.xlsx", sheet = 1)
start_date <- as.Date("1961-01-01")
date_index <- seq(start_date, by = "day", length.out = nrow(catchment_data_sheet2))
catchment_zoo <- zoo(catchment_data_sheet2, order.by = date_index)

# Get year from time index
year <- as.numeric(format(index(catchment_zoo), "%Y"))

# --------- 🔄 MODEL SETUP: USE gr4j_roryVar INSTEAD OF gr4j ----------
mod <- hydromad(
  catchment_zoo,
  sma = "gr4j_roryVar",  # <-- NEW SMA FUNCTION
  routing = "gr4jrouting",
  x1 = c(10, 1000), x2 = c(-5, 5), x3 = c(10, 300), x4 = c(1, 4),
  x5 = c(0.01, 1), x6 = c(0.01, 1), x7 = c(0.01, 1)
)

# Set DE optimization control
Nparams <- 7
paramMult <- 10
ctrl <- DEoptim::DEoptim.control(NP = Nparams * paramMult)

# --------- Objective Function with Weighted Dry & Non-Dry NSEs ----------
F = hydromad.stats(
  'viney' = function(Q, X, dry_years, beta, ...) {
    Q_dry <- Q[(year %in% dry_years), drop = FALSE]
    X_dry <- X[(year %in% dry_years), drop = FALSE]
    Q_nondry <- Q[!(year %in% dry_years), drop = FALSE]
    X_nondry <- X[!(year %in% dry_years), drop = FALSE]
    
    Monthly_NSE_dry <- hmadstat('r.sq.monthly')(Q_dry, X_dry)
    Monthly_NSE_nondry <- hmadstat('r.sq.monthly')(Q_nondry, X_nondry)
    
    NSE_combined <- beta * Monthly_NSE_dry + (1 - beta) * Monthly_NSE_nondry
    return(NSE_combined)
  }
)$viney

# --------- LOOP THROUGH BETA VALUES ----------
beta_values <- seq(0, 1, by = 0.1)
fits_list <- list()
pareto_results <- data.frame(beta = numeric(), nse_dry = numeric(), nse_nondry = numeric(), nse_combined = numeric())

for (beta in beta_values) {
  cat("Running fit for beta =", beta, "\n")
  
  fit <- fitByDE(mod, objective = ~ F(Q, X, dry_years = dry_years, beta = beta))
  
  fitted_Q <- fitted(fit)
  observed_Q <- observed(fit)
  
  fitted_dry <- fitted_Q[year %in% dry_years]
  obs_dry <- observed_Q[year %in% dry_years]
  
  fitted_nondry <- fitted_Q[!(year %in% dry_years)]
  obs_nondry <- observed_Q[!(year %in% dry_years)]
  
  nse_dry <- hmadstat('r.sq.monthly')(obs_dry, fitted_dry)
  nse_nondry <- hmadstat('r.sq.monthly')(obs_nondry, fitted_nondry)
  nse_combined <- beta * nse_dry + (1 - beta) * nse_nondry
  
  pareto_results <- rbind(pareto_results, data.frame(
    beta = beta,
    nse_dry = nse_dry,
    nse_nondry = nse_nondry,
    nse_combined = nse_combined
  ))
  
  fits_list[[paste0("beta_", beta)]] <- fit
}

# --------- PLOTTING ---------
pareto_sorted <- pareto_results[order(pareto_results$nse_dry), ]

plot(pareto_sorted$nse_dry, pareto_sorted$nse_nondry,
     type = "b", pch = 16,
     xlab = "Dry Period NSE",
     ylab = "Non-Dry Period NSE",
     main = "Pareto Curve — Dry vs Non-Dry Fit (GR4J Rory)")

text(pareto_sorted$nse_dry, pareto_sorted$nse_nondry,
     labels = round(pareto_sorted$beta, 2),
     pos = 4, cex = 0.7)


