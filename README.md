# ERIS: Polarization, Antagonism and Community Boundaries

Implementation related to the following article:
Guyot, A., Gillet, A., Leclercq, É., & Cullot, N. (2022, May). ERIS: an approach based on community boundaries to assess polarization in online social networks. In International Conference on Research Challenges in Information Science (pp. 88-104). Cham: Springer International Publishing. Link: https://hal.science/hal-03889719/file/RCIS_2022_Polarization_FINAL.pdf.

Detection and characterization of polarization are of major interest in Social Network Analysis, especially to identify conflictual topics that animate the interactions between users. As gatekeepers of their community, users in the boundaries significantly contribute to its polarization. We propose ERIS, a formal graph approach relying on community boundaries and users’ interactions to compute two metrics: the community antagonism and the porosity of boundaries. These values assess the degree of opposition between communities and their aversion to external exposure, allowing an understanding of the overall polarization through the behaviors of the different communities. We also propose an implementation based on matrix computations, freely available online. Our experiments show a significant improvement in terms of time efficiency in comparison to existing solutions. In the article, we apply our proposal on real data harvested from Twitter with a case study about the vaccines and the COVID-19.

You can find extra resources related to execution times on large graphs in experiment_complexity/.

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
