# This part of the code requires the 'Matrix' library.
if (!require("Matrix")) install.packages("Matrix"); library(Matrix)
require("igraph")

#' Compute the polarization results from an adjacency matrix and a community matrix.
#' 
#' @param Ma An adjacency matrix.
#' @param Mc A community matrix (each row = vertex, each column = community, cell = 1 if vertex in community or 0 otherwise).
#' @param vertices_names A string vector containing the name of each vertex (in the same order as Ma and Mc's rows).
#' @param community_names A string vector containing the name of each community (in the same order as Mc's columns).
#' 
#' @return the antagonism matrix.
compute_polarization_metrics <- function(Ma, Mc, vertices_names = NULL, community_names = NULL, debug = FALSE) {
  nb_comm = ncol(Mc)    # Community count
  
  # Compute the results.
  NMc = 1*!Mc    # The inverse of the binary community matrix.
  McT = t(Mc)    # The transpose of the binary community matrix.
  
  Md = Ma %*% Mc    # Matrix containing the sums of degrees toward each community for each vertex.
  
  I = 1 * (Md == 0)    # A mask of the previous matrix.
  NI = (1 * !I)        # The inverse of this mask.
  
  Ii = vector(mode = "list", length = nb_comm); Mdsi = vector(mode = "list", length = nb_comm); mMdsi = vector(mode = "list", length = nb_comm); Mvani = vector(mode = "list", length = nb_comm); Mani = vector(mode = "list", length = nb_comm); 
  for(i in 1:nb_comm) {
      # A mask revealing the positions of the internal vertices for the community.
      Ii[[i]] = (I * NMc) * Mc[,i]
      # Matrix containing the sums of degrees toward each internal vertex of the community.
      Mdsi[[i]] = (Ma %*% Ii[[i]]) * Mc[,i] * NI
      # Mask of the previous matrix.
      mMdsi[[i]] = 1 * (Mdsi[[i]] != 0)
      # Matrix containing the antagonism scores toward each community for each vertex of the community.
      Mvani[[i]] = (Mdsi[[i]]/(Mdsi[[i]] + Md)) - 0.5; Mvani[[i]] = Mvani[[i]] * (Mdsi[[i]] > 0); Mvani[[i]][is.na(Mvani[[i]])] = 0
      # Vector containing the antagonism scores toward each community for the current community.
      Mani[[i]] = (McT[i,] %*% Mvani[[i]]) / (McT[i,] %*% mMdsi[[i]]); Mani[[i]][is.na(Mani[[i]])] = 0
  }
  
  # Gather the vectors into matrices.
  Man = do.call(rbind, Mani);
  
  # Build the result object.
  res = list(antagonism_matrix = Man)
  
  return (res)
}

#' Compute the polarization results from a graph object created with the igraph library.
#' The weights must be contained in an edge attribute called 'weight'.
#' 
#' @param graph A graph.
#' @param community_attr The vertex attribute containing its community(ies).
#' @param vnames_attr The vertex attribute containing its name/identifier.
#' 
#' @return the antagonism matrix.
compute_polarization_metrics_graph <- function(graph, community_attr = "community", vnames_attr = "name") {
  # Build the adjacency matrix.
  Ma = as_adjacency_matrix(graph, attr = "weight")
  
  # Build the community matrix.
  community_attributes = get.vertex.attribute(graph, community_attr)
  community_names = unique(community_attributes)
  vertices_names = get.vertex.attribute(graph, vnames_attr)
  Mc = t(sapply(community_attributes, function(x) as.numeric(x == community_names)))
  colnames(Mc) = community_names; rownames(Mc) = vertices_names
  
  # Compute and return the polarization results.
  return (list(Mc = Mc, Ma = Ma))
}