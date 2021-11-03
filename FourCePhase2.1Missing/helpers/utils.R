concat <- function(x, y){
  paste0(x, ' (', round(y, 3)*100, '%)')
}

concat_median <- function(med, mi, ma){
  paste0(med, ' [', mi, ', ', ma, ']')
}

concat_mean <- function(mea, s, acc = 0){
  paste0(round(mea, acc), ' (', round(s, acc), ')')
}

right_join0 <- function(x, y, fill = 0L, ...){
  z <- right_join(x, y, ...)
  tmp <- setdiff(names(z), names(y))
  tidyr::replace_na(z, setNames(as.list(rep(fill, length(tmp))), tmp))
}

install_required <- function(required_pkgs){
  new_packages <- required_pkgs[!(required_pkgs %in% installed.packages()[,"Package"])]
  
  if(length(new_packages)) 
    install.packages(new_packages)
}
