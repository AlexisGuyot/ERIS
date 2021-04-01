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
  #if(!is.list(community_members)) community_membership = as.list(community_members)
  #else community_membership = community_members
  structural_matrix = matrix(list(0), nrow = vertex_count, ncol = community_count)
  
  # Detect Internals
  for(v in 1:vertex_count) {
    communities_v = community_membership[[v]]
    for(cv in 1:length(communities_v)) {
      community_v = communities_v[[cv]]
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
  
  # Detect Boundaries
  for(v in 1:vertex_count) {
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
            for(co in 1:length(structural_matrix[neighbor,])) if(structural_matrix[neighbor,][[co]][[cn]] == 0) internal_neighbors = c(internal_neighbors, co)
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
  
  antagonism_matrix <- matrix(0, nrow = community_count, ncol = community_count)
  boundaries_count <- matrix(0, nrow = community_count, ncol = community_count)
  boundaries <- data.frame()
  internals <- data.frame()
  
  formula <- NULL
  if(is.null(adjacency_matrix)) formula <- formula_for_unweighted
  else formula <- formula_for_weighted
  
  for(v in 1:nrow(structural_matrix)) {
    communities_i = community_membership[[v]]
    for(ci in 1:length(communities_i)) {
      community_i = communities_i[[ci]]
      
      communities_j = numeric()
      vertices_internals = numeric()
      for(c in 1:length(structural_matrix[v,])) 
        if(structural_matrix[[v,c]][[ci]] == 2) communities_j = c(communities_j, c)
        else if(structural_matrix[[v,c]][[ci]] == 0) vertices_internals = c(vertices_internals, c)
          
      if(length(communities_j) > 0) for(community_j in communities_j) {
        Ebv = 0
        Eiv = 0
        for(neighbor in adjacency_list[[v]])
          if(community_j %in% community_membership[[neighbor]]) Ebv = formula(Ebv, adjacency_matrix, v, neighbor)
          else {
            index = which(community_membership[[neighbor]] == community_i)
            if(length(index) > 0 && structural_matrix[neighbor,community_j][[index]] == 0) Eiv = formula(Eiv, adjacency_matrix, v, neighbor)
          }
        
        antagonism_matrix[[community_i, community_j]] = antagonism_matrix[[community_i, community_j]] + (Eiv/(Eiv+Ebv) - 0.5)
        boundaries_count[[community_i, community_j]] = boundaries_count[[community_i, community_j]] + 1
        boundaries = rbind(boundaries, data.frame(vertex = v, degree = length(adjacency_list[[v]]), community_vertex = communities_names[community_i], other_community = communities_names[community_j], Pv = (Eiv/(Eiv+Ebv) - 0.5)))
      }
      
      if(length(vertices_internals) > 0) for(intern in vertices_internals) internals = rbind(internals, data.frame(vertex = v, degree = length(adjacency_list[[v]]), community_vertex = communities_names[community_i], other_community = communities_names[intern]))
    }
  }
  
  return (list(boundaries = boundaries, internals = internals, antagonism_matrix = (antagonism_matrix / ifelse(boundaries_count==0, 1, boundaries_count))))
}

#' Calculate the porosity of the boundaries
#' 
#' @param boundaries Boundaries extracted by the function 'build_antagonism_matrix'
#' @param community_membership A list with for each vertex of the graph its community's index
#' 
#' @return A data frame with for each community its porosity value
porosity = function(boundaries, community_membership) {
  percent = function(number_double) { return (paste(round(number_double * 100), "%", sep=''))}
  communities = names(sort(table(unlist(community_membership)), decreasing = TRUE))
  res = data.frame()
  
  for(communauty in communities) {
    n = p = 0
    lines = boundaries[which(boundaries$community_vertex == communauty),]
    if(nrow(lines) > 0) for(i in 1:nrow(lines)) {
      line = lines[i,]
      if(line$Pv <= 0) n = n+1
      else p = p+1
    }
    res = rbind(res, data.frame("community" = communauty, "porosity" = percent(as.double(n/(n+p)))))
  }
  
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
  
  apply_temporary_names <- function(x) { 
    tmp = integer()
    for(xx in 1:length(x)) tmp[[xx]] = as.integer(which(communities_names == x[[xx]]))
    print(tmp)
  }
  c_membership = lapply(community_membership, apply_temporary_names)
  print(c_membership)

  # Build structural list
  structural_matrix <- build_structural_matrix(community_count, vertex_count, adjacency_list, c_membership)
  print(structural_matrix)
  
  # Build antagonism matrix
  res_object = build_antagonism_matrix(structural_matrix, adjacency_list, adjacency_matrix, c_membership, community_count, communities_names)
  colnames(res_object$antagonism_matrix) = communities_names
  rownames(res_object$antagonism_matrix) = communities_names
  print(res_object)
  
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
