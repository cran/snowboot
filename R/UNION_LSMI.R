#' Snowball sampling with multiple inclusion.
#'
#' The function will creat a list of LSMI objects. The function is primairly
#' used in cross-validation.
#' @note see \code{\link{LSMI}.}
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @param n.seeds A numeric vector of seeds for snowball sampling.
#'    It must be a positive integers.
#' @param n.neigh A numeric vector of waves to be sampled around each seed in LSMI.
#'    For example, n.neigh = 0 corresponds to seeds only, and n.neigh = 1
#'    corresponds to sampling seeds and their first neighbors).
#'    Note that the algorithm allows for multiple inclusions.
#' @inheritParams Oempdegreedistrib
#'
#' @return A list containing LSMI objects
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- Union_LSMI(net, n.seeds = c(5,10), n.neigh = 1:2)
Union_LSMI=function(net,n.seeds,n.neigh,seeds=NULL){
  max_seeds<-n.seeds[which.max(n.seeds)]
  sequence_seeds<-sort(n.seeds,decreasing = T)
  num_seeds<-length(n.seeds)
  num_waves<-length(n.neigh)
  # Seed selection: is without replacement and at random
  if (is.null(seeds)) {
    seed0 <- sort(sample(1:length(net$degree), max_seeds, replace = FALSE))
    seeds<-t(as.matrix(seed0))
  } else {
    #in this part,we must make sure that the input of seeds must be the max-num of seeds
    seed0 <- as.vector(seeds)
  }
  # do n.seeds which is the biggest first
  Union_LSMI_output<-as.list(rep(NA,length(n.seeds)*length(n.neigh)))
  for(i in 1:num_waves){
  Union_LSMI_output[[i]]<-LSMI(net,n.seeds = max_seeds,n.neigh = n.neigh[i],seeds = seeds)
  }
  # do n.seeds beside the biggest one and make sure that the length of the n.seeds larger than one
  if(length(n.seeds)>1){
  seed1<-seed0
  for(j in 2:num_seeds){
    #sampling from previous n.seeds-id
    seed1<-sample(seed1,sequence_seeds[j],replace = FALSE)
    seed1<-t(as.matrix(seed1))
    for(k in 1:num_waves){
      Union_LSMI_output[[k+(j-1)*num_waves]]<-LSMI(net,n.seeds = length(seed1),n.neigh = n.neigh[k],seeds = seed1)
    }
  }
  }else{
  Union_LSMI_output=Union_LSMI_output
  }
  #
  return(Union_LSMI_output)

}
