# Polarization, Antagonism and Community Boundaries

Guerra's method extension to conclude about polarization thanks to antagonism between communities : Guerra, P., Meira Jr, W., Cardie, C., & Kleinberg, R. (2013, June). A measure of polarization on social media networks based on community boundaries. In Proceedings of the International AAAI Conference on Web and Social Media (Vol. 7, No. 1).

## To use the code in R

1. Include the script
  source('polarization.R')
  
2. Use the code

2.1 With igraph library :
- Build your graph.
- Use a community detection algorithm.
- Extract communities membership (function membership in igraph)
- Call the function graph_polarization(your_graph, community_membership, weight_attr_name = your_weight_attr_name)
- Use the returned indicatives ($antagonism_matrix, $boundaries, $internals, $porosity)

2.2 Without igraph library :
- Build your graph's adjacency list.
- Use a community detection algorithm.
- Extract community membership for each vertex.
- Call the function polarization(your_adjacency_list, community_membership)
- Use the returned indicatives ($antagonism_matrix, $boundaries, $internals, $porosity)
