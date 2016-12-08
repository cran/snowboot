###################################################################
#modify some places in old version
B.EmpDistrib <-function(net, n.seeds, n.neigh, sam.size=1, n.boot,
                        method = "w", seeds0=NULL,otherNetParameters=FALSE, grow = T){
  #sam.size (==1 for LSMI1) is the number of different samples taken from the network for each i and j
  #otherNetParameters is true if intervals and fallins for the rest of the parmeters
  #  (other than mean) are required.
  Obs.distrib.out <- empd <- as.list(rep(NA, length(n.seeds)*length(n.neigh)))
  #insert the Union_LSMI function to set the unchanged n.seeds id
  union_seeds<-Union_LSMI(net,n.seeds,n.neigh,seeds=seeds0)
  #num of n.seeds combination,is a vector rather than scalar
  num_n.seeds=length(n.seeds)
  num_n.neigh=length(n.neigh)
  #the sequence is expansion
  sequence_n.seeds=sort(n.seeds)
  counter <- 1
  for(i in n.seeds){
    for(j in n.neigh){
      if(grow == T){
        seeds<-union_seeds[[(i/i)+num_n.neigh*(num_n.seeds-which(sequence_n.seeds==i))]][[1]]
      } else {
        seeds = NULL
      }
      if(j==0){
        Obs.distrib<-Oempdegreedistrib(net, n.seeds=i, n.neigh=j, num.sam=sam.size,seeds=seeds)
        #Do we need use lower character tmp instead ?
        TMP <- Obs.distrib$seeds1
      }else{
        seeds<-union_seeds[[(i/i)+num_n.neigh*(num_n.seeds-which(sequence_n.seeds==i))]][[1]]
        Obs.distrib<-Oempdegreedistrib(net, n.seeds=i,n.neigh=j, num.sam=sam.size, seeds=seeds)
      }
      Oparam<-OparametersEst(Obs.distrib)
      #B.distrib<-bootdeg(Obs.distrib, num.sam=sam.size,n.boot=n.boot)
      #return(B.distrib)
      #browser()
      tmp <- bootdeg(Obs.distrib, num.sam=sam.size,n.boot=n.boot, method = method)$empd[[1]]
      Obs.distrib.out[[counter]] <-Obs.distrib
      empd[[counter]] <- tmp
      counter <- counter+1
    }
    #return(B.distrib)
  }
  #the first for loop end
  return(list(Obs.distrib.out=Obs.distrib.out, empd = empd))
}
###################################################################
combineLSMINodes <- function(bootEmpD){
  #function combines the unodes(or seed1) from the elements of
  #Obs.distrib.out object, which is inside the bootEmpD list
  nodes <- NULL
  for(i in 1:length(bootEmpD$Obs.distrib.out)){
    if("nodes_of_LSMI"%in%names(bootEmpD$Obs.distrib.out[[i]])){
      tmp=bootEmpD$Obs.distrib.out[[i]]$nodes_of_LSMI
    } else tmp=bootEmpD$Obs.distrib.out[[i]]$seeds1
    nodes=c(nodes,tmp)
  }
  unlist(nodes)
}
########################################cross_validation###########
#' A function that uses cross-validation to select seed-wave combination for
#' estimation of a degree's frequency.
#'
#' The function's inputs are a network, a vector of possible seed sample-sizes,
#' a vector of possible waves, and a few tuning parameters. The output will
#' contain the best seed-wave combination for each degree and the width of the
#' 95 percent bootstrap confidence intervals at each degree for
#' the best seed-wave combination.
#' @note Only one LSMI per seed-wave combination is currently supported.
#' @references Efron, B. (1979). Bootstrap methods: another look at the
#'  jackknife. The annals of Statistics, 1-26.
#' @references Thompson, M. E., Ramirez Ramirez, L. L., Lyubchich, V. and
#'  Gel, Y. R. (2015), Using the bootstrap for statistical inference
#'  on random graphs. Can J Statistics. doi: 10.1002/cjs.11271
#' @param alpha Desided type I error for bootstrap confidence intervals, which
#'  are obtained using the quantile method.
#' @param network A network object that is list containing:
#'  \describe{
#'    \item{edges}{The edgelist of the network. A two column
#'      \code{matrix} where each row is an edge.}
#'    \item{degree}{The degree sequence of the network, which is
#'      an \code{integer} vector of length n.}
#'    \item{n}{The network order.}
#'  }
#'    The object can be created by \code{\link{local.network.MR.new5}} or
#'    it can be imported.
#' @param n.seeds A numeric vector for the different sample sizes of seed to use
#'  in cross-validation.
#' @param n.neigh A numeric vector for the different waves to use
#'  in cross-validation.
#' @param n.boot The number of bootstrap sample.
#' @param method Can be either "w" for weighted bootstrap or "nw" for
#'    non-weighted bootstrap. "w" is recommended and set as the default method.
#' @param proxyRep The number of time to sample a proxy. Default is 19.
#' @param proxyOrder The size of the proxy sample. Default is 30.
#' @return A list consisting of
#'  \item{selected_seed_wave}{A matrices that provides
#'    the best seed-wave combinations (obtained via cross-validation) for
#'    the respective estimation method.}
#'  \item{selected_seed_wave}{A vector of length 2 that provides
#'    the bootstrap confidence intervals for the estimated mean degree
#'    using the best seed-wave combinations (see above).}
#' @export
#' @examples
#' net <- artificial_networks[[1]]
#' a <- cross_validation_mean(network = net, n.seeds = c(10, 20, 30), n.neigh = c(1, 2),
#'  n.boot = 200, method = "w")

cross_validation_mean <- function(network, n.seeds, n.neigh, n.boot,
                                  method = "w", alpha = .05,proxyRep=19,proxyOrder=30){
  sam.size = 1
  n.seeds <- sort(n.seeds)
  n.neigh <- sort(n.neigh)
  net_order <- network$n
  #make bootEmpD list for seed-wave combos
  bootEmpD <- B.EmpDistrib(net = network, n.seeds = n.seeds, n.neigh = n.neigh,
                           sam.size = sam.size, n.boot = n.boot, method = method)
  fallin.proxy <- array(0, c(length(n.seeds)*length(n.neigh), proxyRep))
  used <- unique(combineLSMINodes(bootEmpD))
  count <- 1
  res<-matrix(NA,length(n.seeds)*length(n.neigh),2)
  for(i in 1:length(n.seeds)){
    # i=1
    for(j in 1:length(n.neigh)){
      # j=1
      # build proxy from bootEmpD$Obs.empd.out
      tmp <- bootEmpD$empd[[count]][[1]]
      values <- bootEmpD$Obs.distrib.out[[count]]$values[[1]]
      est_means <- rowSums(tmp*rep(values,each=n.boot))
      bootCI_mean <- stats::quantile(est_means, c((alpha/2), 1-(alpha/2)))
      res[(i-1)*length(n.neigh)+j,]<-bootCI_mean
      count <- count+1
    }
  }
#calculate the fallin proxy#
for(p in 1:(length(n.seeds)*length(n.neigh))){
  for(k in 1:proxyRep){
  # k=1
  proxyNodes <- sample(used, proxyOrder, replace = F)
  proxy_mean <- mean(network$degree[proxyNodes])
  fallin.proxy[p, k] <- ((res[p,1]<proxy_mean) & (proxy_mean<res[p,2]))
}
}
  coverage.proxy <- apply(fallin.proxy, 1, mean, na.rm=T)
#tranfer to matrix and make a transpose which will be convenient to pick up corresponding index
  matrix_coverage.proxy<-t(matrix(coverage.proxy,length(n.neigh),length(n.seeds)))
  colnames(matrix_coverage.proxy)<-n.neigh
  rownames(matrix_coverage.proxy)<-n.seeds
#max_coverage_lable is the result of combinations which shows their column number and row number
  matrix_coverage.proxy<-abs(matrix_coverage.proxy-0.95)
  max_coverage_label<-which(matrix_coverage.proxy == min(matrix_coverage.proxy), arr.ind = TRUE)
#then judge whether only one seed-wave combination has max coverage.proxy#check what if the cloest is a tie
if(dim(max_coverage_label)[1]==1){
  best_coverage_combination<-data.frame(n.seeds.num=n.seeds[max_coverage_label[1,1]],n.neigh.num=n.neigh[max_coverage_label[1,2]])
}else{
  dimnames(max_coverage_label)<-NULL
  #####when all combination have the same number of seeds
  bigger_n.seeds_index<-which(max_coverage_label[,1]!=min(max_coverage_label[,1]), arr.ind = TRUE)
if(length(bigger_n.seeds_index)==0){
  smaller_n.neigh_index<-which(max_coverage_label[,2]==min(max_coverage_label[,2]))
  #####
  best_coverage_combination<-data.frame(n.seeds.num=n.seeds[max_coverage_label[smaller_n.neigh_index,1]],n.neigh.num=n.neigh[max_coverage_label[smaller_n.neigh_index,2]])
}
  #####when only one combination has a relative small number of seeds
  #####add, bigger_n.seeds_index should be put there, because if all combination with coverage which are the most closest to 95% have the same number of seeds
  #####so if put this part above the second if, the result\bigger_n.seeds_index will be an empty set
  #####because insert one column\the definition of bigger_n.seeds_index in this part,I am not enable to use "else if" below.Becasue it should be
  #####empty between if and else\else if.
else if(length(bigger_n.seeds_index)==(length(max_coverage_label[,1])-1)){
  #####over
  coverage_matrix<-max_coverage_label[-c(bigger_n.seeds_index),]
  best_coverage_combination<-data.frame(n.seeds.num=n.seeds[coverage_matrix[1]],n.neigh.num=n.neigh[coverage_matrix[2]])
}else{
  #####
  coverage_matrix<-max_coverage_label[-c(bigger_n.seeds_index),]
  smaller_n.neigh_index<-which(coverage_matrix[,2]==min(coverage_matrix[,2]))
  best_coverage_combination<-data.frame(n.seeds.num=n.seeds[coverage_matrix[smaller_n.neigh_index,1]],n.neigh.num=n.neigh[coverage_matrix[smaller_n.neigh_index,2]])
}
}
  best_bootCI_mean<-as.matrix(res[(((which(n.seeds==best_coverage_combination[[1]]))-1)*length(n.neigh)+(which(n.neigh==best_coverage_combination[[2]]))),])
  rownames(best_bootCI_mean)<-c("2.5%","97.5%")
  best_coverage_value<-coverage.proxy[(((which(n.seeds==best_coverage_combination[[1]]))-1)*length(n.neigh)+(which(n.neigh==best_coverage_combination[[2]])))]
  return(list(best_combination=best_coverage_combination,best_bootCI_mean=best_bootCI_mean,best_proxy_coverage_value=best_coverage_value))
  }
########################################over###########################################

