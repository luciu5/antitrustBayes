## Create the bundled pnb dataset for antitrustBayes examples
## Source: Philadelphia National Bank merger case (1960)
## Single market, single year, 39 banks

load("d:/Projects/mergerBayes/data/supreme_data.RData")

pnb_raw <- simdata[simdata$event_id == 1, ]
yrs <- as.numeric(as.character(pnb_raw$year))
pnb_raw <- pnb_raw[yrs == min(yrs), ]
pnb_raw <- pnb_raw[pnb_raw$margin >= 0.1 & is.finite(1/pnb_raw$margin) & pnb_raw$margin > 0, ]

# Build the example data frame
# quantities = total_deposits (raw, in bank-level units)
# margins = Lerner index (p-c)/p computed from absolute margins
pnb <- data.frame(
  bank       = pnb_raw$bank_name,
  group      = as.character(pnb_raw$tophold),
  prices     = pnb_raw$rate_deposits,
  quantities = pnb_raw$total_deposits,
  margins    = pnb_raw$margin / pnb_raw$rate_deposits,
  stringsAsFactors = FALSE
)

# Sort by quantity (largest first)
pnb <- pnb[order(-pnb$quantities), ]
rownames(pnb) <- NULL

# Acquirer: Philadelphia NB (tophold 298816)
# Target:   Girard Trust (tophold 889018)
attr(pnb, "acquirer") <- "298816"
attr(pnb, "target")   <- "889018"
attr(pnb, "insideSize") <- sum(pnb_raw$total_deposits) * 1e6

usethis::use_data(pnb, overwrite = TRUE)
