args <- commandArgs(TRUE)
stopifnot(length(args)==2L)
source(args[[2L]],local=TRUE)
expr <- parse(args[[1L]])
extract <- function(name) {
 hits <- Filter(function(e) is.call(e) && length(e)>=3L && identical(e[[1L]],as.name('<-')) && identical(e[[2L]],as.name(name)),expr)
 stopifnot(length(hits)==1L);eval(hits[[1L]],envir=environment(extract))
}
for(n in c('ledger_now','ledger_safe_name','ledger_write_row','ledger_start_row','ledger_complete_row','run_job'))extract(n)
ledger_dir <- tempfile('ledger-');dir.create(ledger_dir)
ledger_run_id <- 'root-control'
configs <- data.frame(config_id=c('A/B','A?B'),demand='logit',conduct='bertrand',seed=1103L,layout='single')
jobs <- data.frame(config_index=1:2,structure='stackelberg',job_index=1:2,job_id=paste(configs$config_id,'stackelberg',sep='__'))
extract('planned_jobs')
atomic_save(planned_jobs,file.path(ledger_dir,'planned_jobs.rds'))
panels <- list(list(),list());fit_dir <- tempfile();rerun_errors <- rerun_diagnostics <- FALSE
fit_one <- function(panel,cc,structure,fit_dir,...) {
 paths <- list.files(ledger_dir,pattern='^[0-9].*\\.rds$',full.names=TRUE)
 rows <- lapply(paths,readRDS)
 stopifnot(any(vapply(rows,function(x) identical(x$state,'started') && identical(x$config_id,cc$config_id),logical(1))))
 list(status='diagnostic_failure',error=NULL)
}
r <- run_job(1L)
stopifnot(identical(r$status,'diagnostic_failure'))
row1 <- readRDS(list.files(ledger_dir,pattern='^0001.*\\.rds$',full.names=TRUE))
stopifnot(identical(row1$state,'completed'),identical(row1$result_status,'diagnostic_failure'),is.finite(row1$pid),!is.na(row1$started_at),!is.na(row1$completed_at))
fit_one <- function(...)stop('controlled error')
e <- tryCatch(run_job(2L),error=identity)
stopifnot(inherits(e,'error'))
row2 <- readRDS(list.files(ledger_dir,pattern='^0002.*\\.rds$',full.names=TRUE))
stopifnot(identical(row2$state,'runner_error'),identical(row2$error,'controlled error'),!is.na(row2$completed_at),!identical(row1$id,row2$id))
stopifnot(all(readRDS(file.path(ledger_dir,'planned_jobs.rds'))$state=='not_started'))
cat('Planned/start/completed/error states, start-before-fit, timestamps, PID, and sanitized-name collision isolation PASS; no fits\n')
