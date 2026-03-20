run_cpp_metabolomics_engine <-
function(dataA,
                                        rt_window = 10,
                                        alpha = 0.7,graphmethod="knn",min_cluster_size=10) {

  library(data.table)
  library(igraph)


  setDT(dataA)
  setorder(dataA, time)

  X <- as.matrix(dataA[, -(1:2)])
  X <- scale(X, center = TRUE, scale = FALSE)

  if(graphmethod=="sparsegraph"){
  #calls the C++ sparse graph builder; alpha is the correlation threshold
  res <- build_sparse_graph_parallel(
    X = X,
    rt = dataA$time,
    rt_window = rt_window,
    alpha = alpha
  )


  edges <- data.table(
    from = as.character(res$from),
    to = as.character(res$to),
    weight = res$weight
  )
  }else{

    print("Running nn2...")
    # Replace non-finite values
    X[is.na(X)] <- 0

    nn <- nn2(X, k = min_cluster_size)

    edges <- cbind(
      rep(1:nrow(X), each = 15),
      as.vector(nn$nn.idx)
    )
}
  g <- graph_from_data_frame(
    edges,
    directed = FALSE,
    vertices = data.frame(name = as.character(seq_len(nrow(dataA))))
  )

  g <- simplify(g, edge.attr.comb = "mean")

  avg_deg <- mean(degree(g))
  resolution <- 1 / avg_deg

  ############################################################
  # Leiden clustering
  ############################################################
  Sys.setenv(OMP_NUM_THREADS = 1)

  set.seed(555)
  cl <- cluster_leiden(
    g,
    weights = E(g)$weight,
    resolution = resolution
  )

  dataA[, Module_RTclust := as.character(membership(cl))]

  return(dataA[, .(mz, time, Module_RTclust)])
}
