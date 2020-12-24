require(corpcor)
require(igraph)

critical_r <- function(n, alpha=0.05) {
  #' Calculate the critical correlation coefficient to determine significance
  #'
  #' @param n the number of data observations
  #' @param alpha the level of significance of the t test
  
  df <- n - 2
  critical.t <- qt(alpha/2, df, lower.tail = F)
  critical.r <- sqrt( (critical.t^2) / ( (critical.t^2) + df ) )
  return(critical.r)
}

get_cig <- function (X, alpha) {
  #' Generate a CIG from a time series sample
  #'
  #' @param X time series data
  #' @param alpha the level of significance of the t test
  #' 
  #' @return graph : the CIG 
  #' @return mat : the partial correlation matrix
  #' @return r: the critical correlation coefficient 
  
  # estimate the covariance matrix
  S = cov(X)
  # get significance theshold
  crit.r = critical_r(nrow(X), alpha)
  # build partial correlation matrix
  pcor.mat = corpcor::cor2pcor(S)
  pcor.mat[upper.tri(pcor.mat, diag=TRUE)] <- 0
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
  diag(pcor.mat) <- 1
  # build CIG from partial correlation matrix
  g <- igraph::graph_from_edgelist(el, directed=FALSE)
  return(list("graph"=g, "mat"=pcor.mat, "r"=crit.r))
}

# TODO: add cig2dag function, given time series CIG and list of variables at time t
