# Independent integration of the inherited bounded correlation-prior kernel.
rmax <- .95
rows <- lapply(c(1,2,4),function(eta) {
  raw <- integrate(function(x)
    (1-(rmax*tanh(x))^2)^(eta-1)/cosh(x)^2,
    -40,40,rel.tol=1e-12)$value
  transformed <- integrate(function(r) (1-r^2)^(eta-1),
    -rmax,rmax,rel.tol=1e-12)$value/rmax
  stopifnot(abs(raw-transformed)<1e-11)
  data.frame(eta,raw_integral=raw,transformed_integral=transformed,
    discrepancy=abs(raw-transformed),log_common_evidence_offset=log(raw))
})
tab <- do.call(rbind,rows)
print(tab,digits=12)
stopifnot(abs(tab$raw_integral[tab$eta==2]-2*(1-rmax^2/3))<1e-11)
cat("PASS: coordinate change and analytic eta=2 integral.\n")
