#!/usr/bin/env Rscript

## Stan-free checks for the compact recovery runner. This file never loads
## antitrustBayes, calls bayes_est(), compiles Stan, or creates fit artifacts.

args <- commandArgs(trailingOnly = FALSE)
f <- args[startsWith(args, "--file=")]
if (!length(f)) stop("test must be invoked with Rscript")
runner <- normalizePath(file.path(dirname(sub("^--file=", "", f[[1L]])),
                                  "recovery_runner.R"), mustWork = TRUE)
options(stackelberg.recovery.runner.dir = dirname(runner))
env <- new.env(parent = globalenv())
sys.source(runner, envir = env)

settings <- env$recovery_primary_settings()
stopifnot(identical(settings$chains, 2L), identical(settings$iter, 1000L),
          identical(settings$warmup, 500L),
          identical(settings$control$max_treedepth, 13L),
          identical(settings$control$adapt_delta, .95),
          identical(settings$smoke, FALSE))

plan <- env$recovery_make_plan()
stopifnot(nrow(plan) == 72L, all(plan$truth_structure %in% c("stackelberg", "fringe")),
          identical(plan$case_id[[1L]], "stackelberg__logit__bertrand__simultaneous__1103__single"),
          all(plan$layout[plan$seed == 1103L] == "single"),
          all(plan$layout[plan$seed == 2207L] == "multi_leader"),
          all(plan$layout[plan$seed == 3301L] == "multi_follower"),
          all(c("stackelberg", "fringe") %in% plan$truth_structure))

cli_default_layout <- env$recovery_parse_cli(c(
  "--demand=logit", "--conduct=bertrand"
))
stopifnot(is.null(cli_default_layout$layout))
cli_plan <- env$recovery_make_plan(
  cli_default_layout$demand, cli_default_layout$conduct,
  cli_default_layout$structure, c(1103L, 2207L, 3301L),
  cli_default_layout$layout,
  truth_structure = cli_default_layout$truth_structure
)
stopifnot(nrow(cli_plan) == 18L,
          identical(unique(cli_plan$layout[cli_plan$seed == 1103L]), "single"),
          identical(unique(cli_plan$layout[cli_plan$seed == 2207L]), "multi_leader"),
          identical(unique(cli_plan$layout[cli_plan$seed == 3301L]), "multi_follower"))

one <- env$recovery_make_plan("logit", "cournot", "fringe", 2207L,
                              "multi_leader", truth_structure = "stackelberg",
                              single_case = TRUE)
stopifnot(nrow(one) == 1L, one$structure == "fringe", one$seed == 2207L)
bad <- tryCatch(env$recovery_make_plan("logit", "cournot", "fringe", 2207L,
                                       c("single", "multi_leader"), single_case = TRUE),
                error = identity)
stopifnot(inherits(bad, "condition"))

fixture_panel <- list(
  data = data.frame(prices = 1:2, quantities = 1:2, margins = c(.2, .3),
                    market = "m", year = c("1", "1"), group = c("A", "B")),
  parameters = list(outside_price = 1.3))
args_exact <- env$recovery_fit_args(fixture_panel,
  list(demand = "logit", conduct = "bertrand", structure = "fringe",
       truth_structure = "stackelberg", seed = 1103L, layout = "single",
       case_id = "x"), settings)
stopifnot(identical(args_exact$partition_prior, "ordered_k"),
          identical(args_exact$k_support, 0:5),
          identical(args_exact$k_prior, rep(1 / 6, 6)),
          identical(args_exact$k_scope, "market"),
          identical(args_exact$stan_options$chains, 2L),
          identical(args_exact$stan_options$iter, 1000L),
          identical(args_exact$stan_options$warmup, 500L),
          identical(args_exact$stan_options$control$max_treedepth, 13L),
          is.integer(args_exact$stan_options$seed), args_exact$stan_options$seed > 0L)
args_stack <- env$recovery_fit_args(fixture_panel,
  list(demand = "ces", conduct = "cournot", structure = "stackelberg",
       truth_structure = "fringe", seed = 2207L, layout = "multi_leader",
       case_id = "y"), settings)
stopifnot(identical(args_stack$k_support, 1:5),
          identical(args_stack$k_prior, rep(1 / 5, 5)))

dgp <- env$recovery_load_dgp()
panel <- dgp("logit", "bertrand", 1103L, "single", normalize_bayes = TRUE,
             partition = "ordered_k", k = 2L, structure = "stackelberg")
env$recovery_validate_dgp(panel)
stopifnot(length(panel$truth) == 8L,
          identical(vapply(panel$truth, function(x) as.integer(x$k_realized), integer(1)),
                    rep(2L, 8L)),
          all(vapply(panel$truth, function(x)
            identical(x$fixed_rank, panel$truth[[1L]]$fixed_rank), logical(1))))
fringe_panel <- dgp("logit", "bertrand", 1103L, "single", normalize_bayes = TRUE,
                    partition = "ordered_k", k = 2L, structure = "fringe")
env$recovery_validate_dgp(fringe_panel, "fringe")

source_text <- paste(readLines(runner, warn = FALSE), collapse = "\n")
stopifnot(grepl("--run", source_text, fixed = TRUE),
          grepl("primary", source_text, fixed = TRUE),
          grepl("never claim recovery from smoke settings", source_text, fixed = TRUE),
          !grepl("rstan::sampling", source_text, fixed = TRUE),
          !grepl("saveRDS\\(.*stanfit.*\\.rds", source_text, perl = TRUE),
          grepl("energy__", source_text, fixed = TRUE),
          grepl("true_k_probability", source_text, fixed = TRUE),
          grepl("margin_mae_core", source_text, fixed = TRUE),
          !grepl("fit <<-", source_text, fixed = TRUE),
          grepl("paste(\"ERROR:\", err)", source_text, fixed = TRUE))

cat("PASS: primary settings, selectors, exact BF/CF/BS/CS supports, and DGP round-trip\n")
cat("PASS: no-fit source checks; HMC/Stan sampling is unreachable without explicit --run\n")
