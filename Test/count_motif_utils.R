##' @noRd
remove_degzero_vs <- function(A) {
  nonzero_ind <- (colSums(A) != 0)
  return(A[nonzero_ind,nonzero_ind])
}

count_motif <- function(A, motif) {
  A <- remove_degzero_vs(A)
  count<-switch(motif, 
                vshape = count_vshape(A),
                triangle = count_triangle(A),
                threestar = count_threestar(A)
  )
  return(count)
}
# counting rooted motifs----
count_rooted_motif <- function(A, motif) {
  count <- switch(motif,
                vshape = count_rooted_vshape(A),
                triangle = count_rooted_triangle(A),
                threestar = count_rooted_threestar(A)
                )
  return(count)
}

count_rooted_triangle <- function(A) {
  return(igraph::count_triangles(igraph::graph_from_adjacency_matrix(A, mode = "undirected")))
}

count_rooted_vshape <- function(A) {
  deg_A <- colSums(A)
  A2 <- A%*%A
  diag(A2) <- 0
  count <- colSums(A2) + (deg_A*(deg_A-1))/2 # (#a-?-?) + (#?-a-?)
  return(count)
}
count_rooted_threestar <- function(A){
  deg_A <- colSums(A)
  n <- length(deg_A)
  
  # #(?- a < ??) + #(a - ?< ??)
  ## #(?- a < ??)
  count <- choose(deg_A, 3)
  ## #(a - ?< ??)
  for(i in 1:n){
    #This only works for 0-1 adjacency matrix
    #The original code used colSums(A^2) and A[,i]^2 to compute choose(deg_A - A[,i],2). 
    count[i] <- count[i] + sum(A[,i] * choose(deg_A-A[,i],2)) 
  }
  return(count)
}

# counting 2v rooted motifs----
count_rooted_2V_motif <- function(A, motif) {
  count<-switch(motif,
                vshape = count_rooted_2V_vshape(A),
                triangle= count_rooted_2V_triangle(A),
                threestar = count_rooted_2V_threestar(A))
  return(count)
}

count_rooted_2V_triangle <- function(A) {
  #i--j (=A_{ij}) & i--?--j (=(A^2)_{ij})
  return(A * (A%*%A))
}
count_rooted_2V_vshape <- function(A) {
  #i--?--j (=(A^2)_{ij}) or i--j--?
  deg_A <- colSums(A)
  A2 <- A%*%A
  diag(A2) <- 0
  count <- A2 + A*pmax(0,outer(deg_A-1, deg_A-1, "+"))
  return(count)
}
count_rooted_2V_threestar <- function(A) {
  #i--j < ? & ?>i--j & i--?<?-j (=(A^2)_{ij})
  deg_A <- colSums(A)
  
  ## i--j < ? & ?>i--j 
  count <- A*outer(choose(deg_A - 1, 2), choose(deg_A - 1, 2), "+")
  ##& i--?<?-j (=(A^2)_{ij})
  count <- count + A %*% diag(pmax(deg_A - 2,0)) %*% A
  
  diag(count) <- 0
  
  return(count)
}
# counting motifs----
count_triangle <- function(A) {
  return(sum(igraph::count_triangles(igraph::graph_from_adjacency_matrix(A, mode = "undirected")))/3)
}

count_vshape <- function(A) {
  deg_A <- colSums(A)
  return(sum(choose(deg_A,2)))
}
count_threestar <- function(A) {
  deg_A <- colSums(A)
  return(sum(choose(deg_A,3)))
}

max_count_motif <-function(num_node, motif){
  return(switch(motif,
                triangle = choose(num_node, 3),
                vshape = 3*choose(num_node, 3),
                threestar = 4 * choose(num_node, 4)
  ))
}
max_root_count_motif <-function(num_node, motif){
  return(switch(motif,
                triangle = choose(num_node-1, 2),
                vshape = 3*choose(num_node-1, 2),
                threestar = 4 * choose(num_node-1, 3)
  ))
}
max_root_2V_count_motif <-function(num_node, motif){
  return(switch(motif,
                triangle = num_node-2,
                vshape = 3*pmax(num_node-2, 0),
                threestar = 6 * choose(num_node-2, 2)
  ))
}