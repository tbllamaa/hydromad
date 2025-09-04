# Trying to get GR4J rainfall runoff model working


# load package - not necessary to do each time when working in a project
library(hydromad)

# Load the new gr7jt model. This loads the altered functions into the environment
# source("gr7jt.R")
# hydromad.model("gr7jt", sma = gr7jt.sim)

hydromad.options(pure.R.code = TRUE)  # ensures we use rcode for gr4j, not C

# Load the HydroTestData dataset. Contains P, E and Q values for a timeseries 2000 Jan to March
data(HydroTestData)



# BUILD the model using the hydromad() function; calling the GR4J sma and routing functions
mod0 <- hydromad(
  HydroTestData, 
  sma = "gr7jt", #gr7jt
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

