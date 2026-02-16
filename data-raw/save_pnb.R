load("d:/Projects/mergerBayes/data/supreme_data.RData")

pnb_raw <- simdata[simdata$event_id == 1, ]
yrs <- as.numeric(as.character(pnb_raw$year))
pnb_raw <- pnb_raw[yrs == min(yrs), ]
pnb_raw <- pnb_raw[pnb_raw$margin >= 0.1 & is.finite(1/pnb_raw$margin) & pnb_raw$margin > 0, ]
pnb_raw$shareIn <- pnb_raw$shareIn / sum(pnb_raw$shareIn)

pnb <- data.frame(
  bank    = pnb_raw$bank_name,
  group   = as.character(pnb_raw$tophold),
  prices  = pnb_raw$rate_deposits,
  shares  = pnb_raw$shareIn,
  margins = pnb_raw$margin,
  stringsAsFactors = FALSE
)

pnb <- pnb[order(-pnb$shares), ]
rownames(pnb) <- NULL

save(pnb, file = "d:/Projects/antitrustBayes/data/pnb.rda", compress = "xz")
cat("Saved pnb.rda with", nrow(pnb), "observations\n")
print(head(pnb))
