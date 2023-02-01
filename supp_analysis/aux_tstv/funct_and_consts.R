#' library for compiling
library(compiler)

#' Preset values and functions
##############################
allow_ezero <- TRUE
goodb <- c("a","c","t","g")
badb <- c("n","-")
ts_lit <- list(c("a","g"),c("c","t"))

#' Retain only good bases
#' a base in outgroup
#' segregating w. two bases in ingroup
check_awful <- function(v,outr=out_rows){check <- TRUE
outu <- unique(v[outr])
inu <- unique(v[-outr])
if (length(intersect(inu,goodb))!=2){check <- FALSE} else {
  if (length(intersect(outu,goodb))<1){check <- FALSE}}                              
return(check)}
ccheck_awful <- cmpfun(check_awful)

#' The reporting function
#' Report whether outgroup has no base, 2+ bases, 1base and whether there is -,n  
#' Report whether ingroup has 
#' 2 bases, 1 base or more
#' missing data
#' gaps
count_bases <- function(v,outr=out_rows){
  out <- c(anc="z",base_oc="z",ominus="z",ona="z",
           base_ic="z",iminus="z",ina="z",pol="z",
           nmut="z","iclean"=-1,"ts"="z")
  #v <- as.character(v)
  outc <- table(v[outr])
  inc <- table(v[-outr])
  outu <- intersect(v[outr],goodb)
  inu <- intersect(v[-outr],goodb)
  out["base_oc"] <- as.character(length(outu))
  out["ominus"] <- as.character(unname(outc["-"]))
  out["ona"] <- as.character(unname(outc["n"]))
  if (out["base_oc"]=="1"){out["anc"] <- outu} else {out["anc"] <- NA}
  out["base_ic"] <- as.character(length(inu))
  out["iminus"] <- as.character(unname(inc["-"]))
  out["ina"] <- as.character(unname(inc["n"]))
  out["ts"] <- NA
  if (length(intersect(v[-outr],badb))==0){out["iclean"] <- "1"} else {
    out["iclean"] <- "0"}
  #' Added check for at least 0.5 base info avail for og
  if (!is.na(out["ona"]) | !is.na((out["ominus"]))){
    if (sum(as.numeric(out["ona"]),
            as.numeric(out["ominus"]),na.rm = TRUE)>length(outr)/2){
      out["iclean"] <- "0"}
  }
  if (!is.na(out["anc"]) & out["base_ic"]=="2"){
    if (any(sapply(ts_lit,function(x){setequal(x,inu)}))){
      out["ts"] <- "1"} else {out["ts"] <- "0"}
    ancin <- intersect(out["anc"],inu)
    if (length(ancin)==0){out["pol"] <- "0";out["nmut"] <- NA} else {
      out["pol"] <- "1"
      mutb <- setdiff(inu,out["anc"])
      out["nmut"] <- as.character(unname(inc[mutb]))}} else {
        out["pol"] <- NA; out["nmut"] <- NA}
  return(out)
}
ccount_bases <- cmpfun(count_bases)

#' Make SFS from data table
sfs_as_sfs <- function(tab1){out1 <- rep(-1,nrow(data2[-out_rows,])-1)
for (i in 1:length(out1)){
  out1[i] <- unname(tab1[as.character(i)])}
if (any(is.na(out1))){out1[is.na(out1)] <- 0}
return(out1)}