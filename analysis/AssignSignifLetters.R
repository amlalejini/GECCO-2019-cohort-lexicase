library(ggplot2)
library(igraph)
library(showtext)
library(fmsb)
library(stringr)

AssignSignifLetters <- function(mw) {
  # And now for a bunch of code that figures out to assign significance group labels for the figure
  # Inspired by http://stackoverflow.com/questions/23681709/algorithm-for-automating-pairwise-significance-grouping-labels-in-r
  
  n <- 10
  g <- as.matrix(mw > 0.05)
  g <- cbind(rbind(NA, g), NA)
  g <- replace(g, is.na(g), FALSE)
  g <- g + t(g)
  diag(g) <- 1
  rownames(g) <- 1:n
  colnames(g) <- 1:n
  
  # Load data
  same <- which(g==1)
  topology <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
  topology <- topology[order(topology[[1]]),] # Get rid of loops and ensure right naming of vertices
  g3 <- simplify(graph.data.frame(topology,directed = FALSE))
  
  # Plot graph
  #plot(g3)
  
  # Calcuate the maximal cliques
  res <- maximal.cliques(g3)
  
  # Reorder given the smallest level
  res <- sapply(res, sort)
  res <- res[order(sapply(res,function(x)paste0(sort(x),collapse=".")))]

  ml<-max(sapply(res, length))
  reord<-do.call("order", data.frame(
    do.call(rbind, 
            lapply(res, function(x) c(sort(as.vector(x)), rep.int(0, ml-length(x))))
    )
  ))
  
  res <- res[reord]
  
  lab.txt <- vector(mode="list", n)
  lab <- letters[seq(res)]
  for(i in seq(res)){
    for(j in res[[i]]){
      lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
    }
  }
  lab.txt
}