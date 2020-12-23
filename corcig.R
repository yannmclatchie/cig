require(corpcor)
require(igraph)

critical.r <- function(n, alpha=0.05) {
  #' Calculate the critical correlation coefficient to determin significance
  #'
  #' @param n the number of data observations
  #' @param alpha the level of significance of the t test
  
  df <- n - 2
  critical.t <- qt(1 - (alpha/2), df)
  critical.r <- sqrt((critical.t^2) / ((critical.t^2) + df))
  return(critical.r)
}

cor2cig <- function (S, varnames=NULL, alpha=0.05) {
  #' Generate a CIG from a covariance matrix
  #'
  #' @param S a sample covariance matrix
  #' @param alpha the level of significance of the t test
  
  # build partial correlation matrix
  parcor.mat = corpcor::cor2pcor(S)
  # get significance theshold
  crit.r = critical.r(nrow(S), alpha=alpha)
  # build variable name dictionary
  var_dict <- as.list(c(colnames(S)))
  names(var_dict) <- c(1:nrow(pcor.mat))
  # store CIG edge list
  el = c()
  # check for significant partial correlation
  for (i in 1:nrow(pcor.mat)) {
    for (j in 1:ncol(pcor.mat)) {
      if (abs(pcor.mat[i,j]) > crit.r) {
        el = rbind(el, c(as.character(var_dict[i]), as.character(var_dict[j])))
      }
    }
  }
  # build CIG from partial correlation matrix
  g <- igraph::graph_from_edgelist(el, directed=FALSE)
  return(g)
}
