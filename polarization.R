library(igraph)

#' Build the structural matrix
#' 
#' @param community_count Number of communities in the graph
#' @param vertex_count Number of vertices in the graph
#' @param adjacency_list Graph's adjacency list (with igraph's adj_list format)
#' @param community_membership A list with for each vertex of the graph its community's index (as.list)
#' 
#' @return The structural matrix
build_structural_matrix <- function(community_count, vertex_count, adjacency_list, community_membership) {
  structural_matrix = matrix(list(0), nrow = vertex_count, ncol = community_count)
  
  print("Début Internals")
  
  # Detect Internals
  for(v in 1:vertex_count) {
    if(v %% 1000 == 0) print(sprintf("%s/%s", v, vertex_count))
    communities_v = community_membership[[v]]
    for(cv in 1:length(communities_v)) {
      community_v = communities_v[[cv]]
      if(cv > 1) for(c in 1:community_count) structural_matrix[[v,c]][[cv]] = 0
      structural_matrix[[v, community_v]][[cv]] = 3
      for(neighbor in adjacency_list[[v]]) {
        communities_neighbor = community_membership[[neighbor]]
        for(cn in 1:length(communities_neighbor)) {
          community_neighbor = communities_neighbor[[cn]]
          if(community_v != community_neighbor) structural_matrix[[v, community_neighbor]][[cv]] = 1
        }
      }
    }
  }
  
  print("Début Boundaries")
  
  # Detect Boundaries
  for(v in 1:vertex_count) {
    if(v %% 1000 == 0) print(sprintf("%s/%s", v, vertex_count))
    communities_v = community_membership[[v]]
    for(cv in 1:length(communities_v)) {
      community_v = communities_v[[cv]]
      external_neighbors = vector()
      internal_neighbors = vector()
      for(neighbor in adjacency_list[[v]]) {
        communities_neighbor = community_membership[[neighbor]]
        for(cn in 1:length(communities_neighbor)) {
          community_neighbor = communities_neighbor[[cn]]
          if(community_v != community_neighbor) external_neighbors = c(external_neighbors, community_neighbor)
          else 
            for(co in 1:length(structural_matrix[neighbor,])) if(structural_matrix[[neighbor,co]][[cn]] == 0) internal_neighbors = c(internal_neighbors, co)
        }
      }
      for(community in intersect(external_neighbors, internal_neighbors)) structural_matrix[[v, community]][[cv]] = 2
    }
  }
  
  return (structural_matrix)
}

#' Build the antagonism matrix
#' 
#' @param structural_matrix The structural matrix built by the 'build_structural_matrix' function
#' @param adjacency_list Graph's adjacency list
#' @param adjacency_matrix Graph's adjacency matrix
#' @param community_membership A list with for each vertex of the graph its community's index
#' @param community_count Number of communities in the graph
#' @param communities_names Communities' names
#' 
#' @return A structure containing the antagonism matrix ($antagonism_matrix), the boundaries members ($boundaries) and the internal areas' members ($internals)
build_antagonism_matrix <- function(structural_matrix, adjacency_list, adjacency_matrix, community_membership, community_count, communities_names) {
  formula_for_weighted <- function(x, adjacency_matrix, source, target) { return (x + adjacency_matrix[source, target]) }
  formula_for_unweighted <- function(x, adjacency_matrix, source, target) { return (x+1) }  
  
  # Structural Matrix content analysis
  internals_size = boundaries_size = 0
  int = bound = 1
  for(elt in structural_matrix) for(value in unlist(elt)) 
    if(value == 0) internals_size = internals_size + 1
    else if(value == 2) boundaries_size = boundaries_size + 1

  # Structures initialization
  antagonism_matrix <- matrix(0, nrow = community_count, ncol = community_count)
  boundaries_count <- matrix(0, nrow = community_count, ncol = community_count)
  boundaries <- matrix(nrow = boundaries_size, ncol = 5, dimnames = list(NULL, c("vertex", "degree", "community_vertex", "other_community", "Pv")))
  internals <- matrix(nrow = internals_size, ncol = 4, dimnames = list(NULL, c("vertex", "degree", "community_vertex", "other_community")))

  # Formula choice (with or whithout weight)
  formula <- NULL
  if(is.null(adjacency_matrix)) formula <- formula_for_unweighted
  else formula <- formula_for_weighted
  
  # Antagonism matrix build
  for(v in 1:nrow(structural_matrix)) {
    if(v %% 1000 == 0) print(sprintf("%s/%s", v, nrow(structural_matrix)))
    communities_i = community_membership[[v]]
    for(ci in 1:length(communities_i)) {
      community_i = communities_i[[ci]]
      communities_j = numeric()
      
      # Structural matrix reading
      for(c in 1:length(structural_matrix[v,])) 
        if(structural_matrix[[v,c]][[ci]] == 2) communities_j = c(communities_j, c)
        else if(structural_matrix[[v,c]][[ci]] == 0) { internals[int,] = c(v, length(adjacency_list[[v]]), communities_names[community_i], c); int = int + 1 }

      # Eiv and Ebv calculations
      if(length(communities_j) > 0) for(community_j in communities_j) {
        Ebv = 0
        Eiv = 0
        for(neighbor in adjacency_list[[v]])
          if(community_j %in% community_membership[[neighbor]]) Ebv = formula(Ebv, adjacency_matrix, v, neighbor)
          else {
            index = match(community_i, community_membership[[neighbor]])
            if(!is.na(index) && structural_matrix[neighbor,community_j][[index]] == 0) Eiv = formula(Eiv, adjacency_matrix, v, neighbor)
          }

        # Antagonism calculation
        antagonism_matrix[[community_i, community_j]] = antagonism_matrix[[community_i, community_j]] + (Eiv/(Eiv+Ebv) - 0.5)
        boundaries_count[[community_i, community_j]] = boundaries_count[[community_i, community_j]] + 1
        boundaries[bound,] = c(v, length(adjacency_list[[v]]), communities_names[community_i], communities_names[community_j], (Eiv/(Eiv+Ebv) - 0.5)); bound = bound + 1
      }
    }
  }
  
  return (list(boundaries = as.data.frame(boundaries), internals = as.data.frame(internals), antagonism_matrix = (antagonism_matrix / ifelse(boundaries_count==0, 1, boundaries_count))))
}

#' Calculate the porosity of the boundaries
#' 
#' @param boundaries Boundaries extracted by the function 'build_antagonism_matrix'
#' @param community_membership A list with for each vertex of the graph its community's index
#' 
#' @return A data frame with for each community its porosity value
porosity = function(boundaries, community_membership) {
  percent = function(number_double) { return (paste(round(number_double * 100), "%", sep=''))}
  distinct = function(dataframe, columns) { return (dataframe[!duplicated(dataframe[,columns]),]) }
  
  # Sizes analysis
  communities_sizes = table(unlist(community_membership))
  boundaries_sizes = table(distinct(boundaries, c('vertex', "community_vertex"))$community_vertex)
  communities = names(communities_sizes)
  res = data.frame()
  
  # Porosity calculation
  for(community in communities) {
    n = p = 0
    lines = boundaries[which(boundaries$community_vertex == community),]
    if(nrow(lines) > 0) for(i in 1:nrow(lines)) {
      line = lines[i,]
      if(as.numeric(as.character(line$Pv)) <= 0) n = n+1
      else p = p+1
    }
    if(n != 0 || p != 0) score = as.double(n/(n+p)) else score = 0
    if(communities_sizes[[community]] > 0 && community <= length(boundaries_sizes)) boundary_size = boundaries_sizes[[community]]/communities_sizes[[community]]
    else boundary_size = 0
    res = rbind(res, data.frame("community" = community, "community_size" = communities_sizes[[community]], "boundary_size" = boundary_size, "porosity" = score))
  }
  
  # Result shaping
  res = res[order(res$porosity),]
  res$boundary_size = sapply(res$boundary_size, percent)
  res$porosity = sapply(res$porosity, percent)
  
  return (res)
}

#' Return indicatives to conclude about polarization on a graph
#' 
#' @param adjacency_list Graph's adjacency list
#' @param community_membership A list with for each vertex of the graph its community's index
#' @param adjacency_matrix Graph's adjacency matrix
#' 
#' @return A structure containing the antagonism matrix ($antagonism_matrix), the boundaries members ($boundaries), the internal areas' members ($internals) and the porosity values ($porosity)
polarization <- function(adjacency_list, community_membership, adjacency_matrix = NULL) {
  vertex_count = length(adjacency_list)
  size_communities = table(unlist(community_membership))
  community_count = length(size_communities)
  communities_names = names(size_communities)
  
  apply_temporary_names <- function(x) {  tmp = integer(); for(xx in 1:length(x)) tmp[[xx]] = as.integer(which(communities_names == x[[xx]])); return (tmp) }
  c_membership = lapply(community_membership, apply_temporary_names)
  
  print("Etape 1")

  # Build structural list
  structural_matrix <- build_structural_matrix(community_count, vertex_count, adjacency_list, c_membership)
  
  print("Etape 2")

  # Build antagonism matrix
  res_object = build_antagonism_matrix(structural_matrix, adjacency_list, adjacency_matrix, c_membership, community_count, communities_names)
  colnames(res_object$antagonism_matrix) = communities_names
  rownames(res_object$antagonism_matrix) = communities_names

  print("Etape 3")
  
  # Calculate porosity
  res_object$porosity = porosity(res_object$boundaries, community_membership)
  
  return (res_object)
}

#' Return indicatives to conclude about polarization on a graph built with igraph
#' 
#' @param graph The graph build with igraph
#' @param community_membership A list with for each vertex of the graph its community's index
#' @param adjacency_matrix Graph's adjacency matrix
#' @param weight_attr_name Name of the edges' attribute containing the weight of the link (default to 'weight')
#' 
#' @return A structure containing the antagonism matrix ($antagonism_matrix), the boundaries members ($boundaries), the internal areas' members ($internals) and the porosity values ($porosity)
graph_polarization <- function(graph, community_membership = NULL, adjacency_matrix = NULL, weight_attr_name = "weight") {
  adjacency_list = as_adj_list(graph, mode = "out")
  if(is.weighted(graph)) adjacency_matrix = as_adjacency_matrix(graph, attr = weight_attr_name)
  if(is.null(community_membership)) community_membership = membership(graph)
  
  # Build antagonism matrix and boundaries
  res_object = polarization(adjacency_list, community_membership, adjacency_matrix)
  
  return (res_object)
}
