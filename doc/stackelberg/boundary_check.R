# Extended research checks; source the independent oracle, not production.
source("doc/stackelberg/math_check.R")
rows<-list(); i<-0
for(gamma in c(1.1,1.5,2,6,15)) for(conduct in c("bertrand","cournot")) {
  cs<-cases$multi_follower; cs$gamma<-gamma
  mk<-make_market(cs,"ces",conduct)
  h<-1e-5
  fd<-sapply(seq_along(mk$li),function(j) {
    plus<-minus<-mk$z[mk$li]; plus[j]<-plus[j]+h; minus[j]<-minus[j]-h
    (mk$solve_f(plus)[mk$fi]-mk$solve_f(minus)[mk$fi])/(2*h)
  })
  i<-i+1; rows[[i]]<-data.frame(demand="ces",conduct,gamma,small_share=NA_real_,reaction_fd=max(abs(fd-mk$reaction)),relative_fd=max(abs(fd-mk$reaction))/max(abs(mk$reaction)))
}
for(demand in c("logit","ces")) for(conduct in c("bertrand","cournot")) for(eps in c(1e-2,1e-4,1e-6,1e-8)) {
  cs<-cases$tiny; cs$r<-c(.25,eps,2*eps,3*eps)
  mk<-make_market(cs,demand,conduct)
  central<-function(h) matrix((mk$solve_f(mk$z[mk$li]+h)[mk$fi]-mk$solve_f(mk$z[mk$li]-h)[mk$fi])/(2*h),ncol=1)
  if(conduct=="bertrand"&&eps<=1e-6) {
    d1<-central(.1); d2<-central(.05); d3<-central(.025)
    fd<-(16*(4*d3-d2)/3-(4*d2-d1)/3)/15
  } else fd<-central(1e-5)
  i<-i+1; rows[[i]]<-data.frame(demand,conduct,gamma=if(demand=="ces")4 else NA_real_,small_share=eps,reaction_fd=max(abs(fd-mk$reaction)),relative_fd=max(abs(fd-mk$reaction))/max(abs(mk$reaction)))
}
tab<-do.call(rbind,rows)
write.csv(tab,file.path(out,"boundary_validation.csv"),row.names=FALSE)
print(tab,digits=4)
stopifnot(max(tab$reaction_fd)<2e-7,max(tab$relative_fd)<2e-5)
cat("PASS: gamma-domain and shrinking-follower reaction checks.\n")
