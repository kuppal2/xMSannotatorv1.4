run_cpp_metabolomics_engine <-
function(dataA,
                                        rt_window = 10,
                                        alpha = 0.7,graphmethod="sparse",min_cluster_size=10) {

  library(data.table)
  library(igraph)


  setDT(dataA)
  setorder(dataA, time)

  X <- as.matrix(dataA[, -(1:2)])
  X <- scale(X, center = TRUE, scale = FALSE)

  if(graphmethod=="sparsegraph"){
    print("Running sparse graph clustering...")
  #calls the C++ sparse graph builder; alpha is the correlation threshold
  res <- build_sparse_graph_parallel(
    X = X,
    rt = dataA$time,
    rt_window = rt_window,
    alpha = alpha,
    top_k=min_cluster_size
  )


  edges <- data.table(
    from = as.character(res$from),
    to = as.character(res$to),
    weight = res$weight
  )
  }else{

    print("Running knn graph clustering...")
    # Replace non-finite values
    X[is.na(X)] <- 0
    X[!is.finite(X)] <- 0  # also catch Inf/-Inf

    nn <- nn2(X, k = min_cluster_size)

    edges <- data.frame(
      from   = as.character(rep(1:nrow(X), each = min_cluster_size)),
      to     = as.character(as.vector(nn$nn.idx)),
      weight = 1 / (1 + as.vector(nn$nn.dists))  # convert distance → similarity
    )
}
  g <- graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = data.frame(name = as.character(seq_len(nrow(dataA))))
  )

  g <- simplify(g, edge.attr.comb = "mean")

  avg_deg <- mean(degree(g))
  #resolution <- 1 / avg_deg
  resolution <- 0.05 #1.0 / log1p(avg_deg)
  ############################################################
  # Leiden clustering
  ############################################################
  Sys.setenv(OMP_NUM_THREADS = 1)
  # Remove NaN/NA/Inf edges
  bad_edges <- which(is.nan(E(g)$weight) | is.na(E(g)$weight) | !is.finite(E(g)$weight))
  #cat("Removing", length(bad_edges), "bad edges\n")
  g <- delete_edges(g, bad_edges)

 # E(g)$weight <- pmax(E(g)$weight^2, 0)
  set.seed(555)

   cl <- cluster_infomap(
    g,
    e.weights = E(g)$weight
    #resolution = resolution
  )


  dataA[, Module_RTclust := as.character(membership(cl))]

  return(dataA[, .(mz, time, Module_RTclust)])
}
