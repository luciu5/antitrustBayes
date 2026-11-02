# Independent research oracle for the exact core-fringe (BF/CF) endpoints.
# Compares coordination::core_fringe() at corePre=character() (k=0) and
# corePre=all firms (k=K) directly against unchanged antitrust::Logit /
# LogitCournot / CES / CESCournot objects (the atomic MC/B/C oracle). This is
# the fringe-family counterpart to ordinary_nesting_validation.R, which only
# exercises the Stackelberg leader/follower endpoints.
#
# Run: Rscript doc/stackelberg/core_fringe_nesting_validation.R [path] [out]
args <- commandArgs(TRUE)
path <- if (length(args)) args[1] else "/tmp/stackelberg-coordination"
out <- if (length(args) > 1) args[2] else "doc/stackelberg"
demands <- if (length(args) > 2) args[3] else c("logit", "ces")
# Only coordination/antitrust (pure-R deterministic games) are needed here;
# antitrustBayes itself (Stan) is not exercised by this oracle.
suppressPackageStartupMessages(library(antitrust))
suppressPackageStartupMessages(library(coordination))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
script <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
oracle <- new.env(parent = globalenv())
for (e in parse(file.path(dirname(script), "math_check.R")))
  if (is.call(e) && identical(e[[1]], as.name("<-")) &&
      as.character(e[[2]]) %in% c("cases", "make_market")) eval(e, oracle)

rows <- list()
i <- 0L
for (demand in demands) for (conduct in c("bertrand", "cournot"))
  for (layout in c("single", "multi_leader", "multi_follower")) for (endpoint in c("all", "none")) {
    cc <- oracle$cases[[layout]]
    # Reuse the leaders field as the core-fringe "core" set: endpoint="none"
    # is k=0 (empty core, exact MC); endpoint="all" is k=K (full core, exact
    # B/C). make_market() is only used here to get consistent prices/shares;
    # its own leaders-dependent h/margin construction is not used below.
    core <- if (endpoint == "all") unique(cc$owner) else character()
    mk <- oracle$make_market(cc, demand, conduct, FALSE)
    obj <- coordination::core_fringe(mk$p, cc$r, mk$m, ownerPre = mk$owner,
      corePre = core, demand = demand, conduct = conduct, output = TRUE,
      insideSize = 1000 * sum(cc$r), priceOutside = 1.7,
      alpha = if (demand == "logit") 2.3 else NULL,
      gamma = if (demand == "ces") 4 else NULL)
    ordinary_class <- if (endpoint == "none") {
      # k=0: every firm is fringe/atomistic -> monopolistic competition,
      # not plain Logit/CES Bertrand (which assumes a strategic firm).
      if (demand == "logit") "MonComLogit" else "MonComCES"
    } else if (demand == "logit") {
      if (conduct == "bertrand") "Logit" else "LogitCournot"
    } else if (conduct == "bertrand") "CES" else "CESCournot"
    ordinary <- methods::new(ordinary_class)
    default_control_equ <- ordinary@control.equ
    for (nm in methods::slotNames(ordinary))
      methods::slot(ordinary, nm) <- methods::slot(obj, nm)
    # obj@control.equ only carries price_domain (set by the coordination
    # constructor); restore the base-class solver defaults (tol/maxit/etc.)
    # that the naive slot-copy above overwrote, merging price_domain back in.
    ordinary@control.equ <- utils::modifyList(default_control_equ, obj@control.equ)
    ordinary@priceStart <- obj@prices * 1.02
    suppressWarnings(pp <- antitrust::calcPrices(ordinary))
    mm <- antitrust::calcMargins(ordinary)
    mc <- antitrust::calcMC(ordinary)
    i <- i + 1L
    rows[[i]] <- data.frame(demand, conduct, layout, endpoint,
      atomic_target = if (endpoint == "none") "MC" else if (conduct == "bertrand") "B" else "C",
      margin = max(abs(mm - antitrust::calcMargins(obj))),
      cost = max(abs(mc - obj@mcPre)),
      price = max(abs(pp - obj@pricePre)))
  }
tab <- do.call(rbind, rows)
write.csv(tab, file.path(out, "core_fringe_nesting_validation.csv"), row.names = FALSE)
print(aggregate(cbind(margin, cost, price) ~ demand + conduct + endpoint, tab, max), digits = 5)
stopifnot(max(tab$margin) < 2e-9, max(tab$cost) < 2e-9, max(tab$price) < 2e-5)
cat("PASS: both core-fringe endpoints (k=0 -> MC, k=K -> B/C) agree with ordinary antitrust conduct.\n")
