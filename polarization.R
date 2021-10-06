# This part of the code requires the 'Matrix' library.
if (!require("Matrix")) install.packages("Matrix"); library(Matrix)
require("igraph")

# Path to the R file containing the functions to generate the vizualisation charts.
VizuLoaded = FALSE
if(file.exists("polarization_vizualisation.R")) source("polarization_vizualisation.R")

#' Compute the polarization metrics from an adjacency matrix and a community matrix.
#' 
#' @param Ma An adjacency matrix.
#' @param Mc A community matrix (each row = vertex, each column = community, cell = 1 if vertex in community or 0 otherwise).
#' @param vertices_names A string vector containing the name of each vertex (in the same order as Ma and Mc's rows).
#' @param community_names A string vector containing the name of each community (in the same order as Mc's columns).
#' @param debug A logical scalar, TRUE returns all the transitional matrices, FALSE only the final ones with the results. 
#' 
#' @return if debug == FALSE: the antagonism matrix, the porosity matrix, matrices with the sizes of the boundaries (absolute and relative) the porosity by boundary size chart ; if debug == TRUE, all the matrices computed without fancy names.
compute_polarization_metrics <- function(Ma, Mc, vertices_names = NULL, community_names = NULL, debug = FALSE) {
  nb_comm = ncol(Mc)    # Community count
  csizes = colSums(Mc)  # Community sizes
  
  # Assign auto-generated or not names to vertices and communities.
  if(!is.null(community_names)) colnames(Mc) = community_names 
  else colnames(Mc) = sapply(seq(nb_comm), function(x) paste('C',x, sep=''))
  if(!is.null(vertices_names)) { colnames(Ma) = vertices_names; rownames(Ma) = vertices_names; rownames(Mc) = vertices_names } 
  else { s = sapply(seq(nrow(Mc)), function(x) paste('',x, sep='')); colnames(Ma) = s; rownames(Ma) = s; rownames(Mc) = s }
  
  # Compute the results.
  NMc = 1*!Mc    # The inverse of the binary community matrix.
  McT = t(Mc)    # The transpose of the binary community matrix.
  
  Md = Ma %*% Mc    # Matrix containing the sums of degrees toward each community for each vertex.
  
  I = 1 * (Md == 0)    # A mask of the previous matrix.
  NI = (1 * !I)        # The inverse of this mask.
  
  Ii = vector(mode = "list", length = nb_comm); Mdsi = vector(mode = "list", length = nb_comm); mMdsi = vector(mode = "list", length = nb_comm); Mvani = vector(mode = "list", length = nb_comm); clean_Mvani = vector(mode = "list", length = nb_comm); Mbsi = vector(mode = "list", length = nb_comm); Mbrsi = vector(mode = "list", length = nb_comm); Mani = vector(mode = "list", length = nb_comm); Mpi = vector(mode = "list", length = nb_comm);
  for(i in 1:nb_comm) {
    # A mask revealing the positions of the internal vertices for the community.
    Ii[[i]] = (I * NMc) * Mc[,i]
    # Matrix containing the sums of degrees toward each internal vertex of the community.
    Mdsi[[i]] = (Ma %*% Ii[[i]]) * Mc[,i] * NI
    # Mask of the previous matrix.
    mMdsi[[i]] = 1 * (Mdsi[[i]] != 0)
    # Matrix containing the antagonism scores toward each community for each vertex of the community.
    Mvani[[i]] = (Mdsi[[i]]/(Mdsi[[i]] + Md)) - 0.5; Mvani[[i]] = Mvani[[i]] * (Mdsi[[i]] > 0); Mvani[[i]][is.na(Mvani[[i]])] = 0
    # Same as the last one but only containing the rows corresponding to a boundary vertex of the community.
    clean_Mvani[[i]] = list(Mvani[[i]][apply(mMdsi[[i]], 1, function(x) !all(x==0)),])
    # Vector containing the sizes of all the boundaries of the current community.
    Mbsi[[i]] = colSums(mMdsi[[i]])
    # Vector containing the antagonism scores toward each community for the current community.
    Mani[[i]] = (McT[i,] %*% Mvani[[i]]) / (McT[i,] %*% mMdsi[[i]]); Mani[[i]][is.na(Mani[[i]])] = 0
    # Vector containing the porosity scores toward each community for the current community.
    Mpi[[i]] = ((McT[i,] %*% (1*(Mvani[[i]] < 0))) / replace(Mbsi[[i]], Mbsi[[i]] == 0, 1)) * 100
    # Transform absolute size into relative size.
    Mbrsi[[i]] = Mbsi[[i]] / csizes[[i]] * 100                         
  }
  
  # Gather the vectors into matrices.
  Man = do.call(rbind, Mani); rownames(Man) = community_names
  Mp = do.call(rbind, Mpi); rownames(Mp) = community_names
  Mbs = do.call(rbind, Mbsi); rownames(Mbs) = community_names
  Mbrs = do.call(rbind, Mbrsi); rownames(Mbrs) = community_names
  
  # Create user-friendly vizualisation of the results.
  Vizu = NULL
  if(VizuLoaded) Vizu = compute_polarization_vizualisations(Mc, Man, Mp, Mbs, Mbrs, community_names)
  
  # Build the result object.
  res = list()                                                                                               
  if(!debug) res = list(antagonism_vertices = clean_Mvani, boundary_sizes_matrix = Mbs, boundary_rsizes_matrix = Mbrs, porosity_matrix = Mp, antagonism_matrix = Man, vizualisation = Vizu)
  else res = list(Ma = Ma, Mc = Mc, NMc = NMc, McT = McT, Md = Md, I = I, NI = NI, Ii = Ii, Mdsi = Mdsi, Mvani = Mvani, clean_Mvani = clean_Mvani, Mani = Mani, Mpi = Mpi, Mp = Mp, Mbsi = Mbsi, Mbrsi = Mbrsi, Mbs = Mbs, Mbrs = Mbrs, Man = Man, vizu = Vizu)
  
  return (res)
}

#' Create the charts from polarization results.
#' 
#' @param community_matrix The community matrix (each row = vertex, each column = community, cell = 1 if vertex in community or 0 otherwise).
#' @param antagonism_matrix The antagonism matrix.
#' @param porosity_matrix The porosity matrix.
#' @param boundary_sizes_matrix The boundary sizes matrix.
#' @param boundary_rsizes_matrix The relative boundary sizes matrix.
#' @param community_names A list with exactly one name by community (in the same order as in the matrices).
#' 
#' @return A list containing all the charts generated with ggplot2.
compute_polarization_vizualisations <- function(community_matrix, antagonism_matrix, porosity_matrix, boundary_sizes_matrix, boundary_rsizes_matrix, community_names = NULL) {
  boundary_sizes = `dim<-`(sprintf('%s\n(%s%%)', boundary_sizes_matrix, round(boundary_rsizes_matrix, 1)), dim(boundary_sizes_matrix))
  colnames(boundary_sizes) = colnames(boundary_rsizes_matrix); rownames(boundary_sizes) = rownames(boundary_rsizes_matrix)
  return (list(
    antagonism_matrix = matrix_vizualisation(antagonism_matrix, community_matrix, low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-0.5,0.5), name="Antagonism", decimal = 3, community_names),
    porosity_matrix = matrix_vizualisation(porosity_matrix, community_matrix, low = "white", high = "green4", mid = "darkolivegreen3", midpoint = 50, limit = c(0,100), name = "Porosity (in %)", decimal = 1, community_names),
    boundary_size_matrix = matrix_vizualisation(boundary_rsizes_matrix, community_matrix, low = "white", high = "darkgray", mid = "gray", midpoint = 50, limit = c(0,100), name = "Boundary size", decimal = 2, community_names, boundary_sizes),
    porosity_by_size = porosity_vizualisation2D(porosity_matrix, boundary_sizes_matrix)
  ))
}

#' Compute the polarization metrics from a graph object created with the igraph library.
#' The weights must be contained in an edge attribute called 'weight'.
#' 
#' @param graph A graph.
#' @param community_attr The vertex attribute containing its community(ies).
#' @param vnames_attr The vertex attribute containing its name/identifier.
#' @param debug A logical scalar, TRUE returns all the transitional matrices, FALSE only the final ones with the results.
#' 
#' @return if debug == FALSE: the antagonism matrix, the porosity matrix, the porosity by boundary size chart ; if debug == TRUE, all the matrices computed without fancy names.
compute_polarization_metrics_graph <- function(graph, community_attr = "community", vnames_attr = "name", debug = FALSE) {
  # Build the adjacency matrix.
  if(!is.null(vertex_attr(graph, "weight"))) Ma = as_adjacency_matrix(graph, attr = "weight")
  else Ma = as_adjacency_matrix(graph)
  
  # Build the community matrix.
  community_attributes = get.vertex.attribute(graph, community_attr)
  community_names = unique(community_attributes)
  vertices_names = get.vertex.attribute(graph, vnames_attr)
  Mc = t(sapply(community_attributes, function(x) as.numeric(x == community_names)))
  colnames(Mc) = community_names; rownames(Mc) = vertices_names
  
  # Compute and return the polarization results.
  return (compute_polarization_metrics(Ma, Mc, vertices_names = vertices_names, community_names = community_names, debug = debug))
}