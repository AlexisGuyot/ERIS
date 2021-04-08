library(igraph)

#' Build the structural matrix
#' 
#' @param community_count Number of communities in the graph
#' @param vertex_count Number of vertices in the graph
#' @param adjacency_list Graph's adjacency list
#' @param community_membership A list with for each vertex of the graph its community's index
#' 
#' @return The structural matrix
build_structural_matrix <- function(community_count, vertex_count, adjacency_list, community_membership) {
  structural_matrix = matrix(0, nrow = vertex_count, ncol = community_count)
  
  print("Debut Internals")
  
  # Detect Internals
  for(v in 1:vertex_count) {
    if(v %% 1000 == 0) print(sprintf("%s/%s", v, vertex_count))
    community_v = community_membership[[v]]
    for(neighbor in adjacency_list[[v]]) {
      community_neighbor = community_membership[[neighbor]]
      structural_matrix[[v, community_v]] = 3
      if(community_v != community_neighbor) structural_matrix[[v, community_neighbor]] = 1
    }
  }
  
  print("Debut Boundaries")
  
  # Detect Boundaries
  for(v in 1:vertex_count) {
    if(v %% 1000 == 0) print(sprintf("%s/%s", v, vertex_count))
    community_v = community_membership[[v]]
    external_neighbors = vector()
    internal_neighbors = vector()
    for(neighbor in adjacency_list[[v]]) {
      community_neighbor = community_membership[[neighbor]]
      if(community_v != community_neighbor) external_neighbors = c(external_neighbors, community_neighbor)
      else internal_neighbors = c(internal_neighbors, match(0, structural_matrix[neighbor,]))
    }
    for(community in intersect(external_neighbors, internal_neighbors)) structural_matrix[[v, community]] = 2
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
    if(v %% 1000 == 0) print(sprintf("%s/%s", v, nrow(structural_matrix)))
    community_i = community_membership[[v]]
    communities_j = match(2, structural_matrix[v,])
    if(!is.na(communities_j)) for(community_j in communities_j) {
      Ebv = 0
      Eiv = 0
      for(neighbor in adjacency_list[[v]])
        if(community_membership[[neighbor]] == community_j) Ebv = formula(Ebv, adjacency_matrix, v, neighbor)
      else if(community_membership[[neighbor]] == community_i && structural_matrix[neighbor,community_j] == 0) Eiv = formula(Eiv, adjacency_matrix, v, neighbor)
      
      antagonism_matrix[[community_i, community_j]] = antagonism_matrix[[community_i, community_j]] + (Eiv/(Eiv+Ebv) - 0.5)
      boundaries_count[[community_i, community_j]] = boundaries_count[[community_i, community_j]] + 1
      boundaries = rbind(boundaries, data.frame(vertex = v, degree = length(adjacency_list[[v]]), community_vertex = communities_names[community_i], other_community = communities_names[community_j], Pv = (Eiv/(Eiv+Ebv) - 0.5)))
    }
    
    vertices_internals = match(0, structural_matrix[v,])
    if(!is.na(vertices_internals)) for(intern in vertices_internals) internals = rbind(internals, data.frame(vertex = v, degree = length(adjacency_list[[v]]), community_vertex = communities_names[community_i], other_community = communities_names[intern]))
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
    communities = names(sort(table(community_membership), decreasing = TRUE))
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
  size_communities = table(community_membership)
  community_count = length(size_communities)
  communities_names = names(size_communities)
  
  c_membership = sapply(community_membership, function(x) which(communities_names == x))
  
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
