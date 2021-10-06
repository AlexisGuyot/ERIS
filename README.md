# ERIS: Polarization, Antagonism and Community Boundaries

Description is coming.

## To use the code in R

1. Include the script
  source('polarization.R')
  
2. Use the code

2.1 With igraph library :
- Build your graph.
- Use a community detection algorithm.
- Extract communities membership (function membership in igraph)
- Add to each vertex its community identifier in an attribute (V(g)$community = membership(g))
- Call the function compute_polarization_metrics_graph(your_graph, community_attr = your_community_attr_name, vnames_attr = your_vertex_name_attr)
- Use the returned indicatives ($antagonism_vertices, $boundary_sizes_matrix, $porosity_matrix, $antagonism_matrix, $vizualisation)

2.2 Without igraph library :
- Build your graph's adjacency matrix.
- Use a community detection algorithm.
- Extract community membership for each vertex.
- Build a binary community matrix where 0 means that the vertex heading the row is not in the community heading the column, and 1 otherwise.
- Call the function compute_polarization_metrics(your_adjacency_matrix your_community_matrix, vertices_names = vector_with_all_vertices_names, community_names = vector_with_all_community_names)
- Use the returned indicatives ($antagonism_vertices, $boundary_sizes_matrix, $porosity_matrix, $antagonism_matrix, $vizualisation)
