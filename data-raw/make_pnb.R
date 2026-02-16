## Create the bundled pnb dataset for antitrustBayes examples
## Source: Philadelphia National Bank merger case (1960)
## Single market, single year, 39 banks

load("d:/Projects/mergerBayes/data/supreme_data.RData")

pnb_raw <- simdata[simdata$event_id == 1, ]
yrs <- as.numeric(as.character(pnb_raw$year))
pnb_raw <- pnb_raw[yrs == min(yrs), ]
pnb_raw <- pnb_raw[pnb_raw$margin >= 0.1 & is.finite(1/pnb_raw$margin) & pnb_raw$margin > 0, ]
pnb_raw$shareIn <- pnb_raw$shareIn / sum(pnb_raw$shareIn)

# Build the example data frame
pnb <- data.frame(
  bank    = pnb_raw$bank_name,
  group   = as.character(pnb_raw$tophold),
  prices  = pnb_raw$rate_deposits,
  shares  = pnb_raw$shareIn,
  margins = pnb_raw$margin,
  stringsAsFactors = FALSE
)

# Sort by share (largest first)
pnb <- pnb[order(-pnb$shares), ]
rownames(pnb) <- NULL

# Acquirer: Philadelphia NB (tophold 298816)
# Target:   Girard Trust (tophold 889018)
attr(pnb, "acquirer") <- "298816"
attr(pnb, "target")   <- "889018"
attr(pnb, "insideSize") <- sum(pnb_raw$total_deposits) * 1e6

usethis::use_data(pnb, overwrite = TRUE)
