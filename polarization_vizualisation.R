# This part of the code requires the 'ggplot2' library.
if (!require("ggplot2")) install.packages("ggplot2"); library(ggplot2)

#' Create a heat map to see easier a matrix.
#' 
#' @param matrix_to_display The matrix to display as a heat map.
#' @param Mc A community matrix (each row = vertex, each column = community, cell = 1 if vertex in community or 0 otherwise).
#' @param low A string informing the color for the lowest value.
#' @param high A string informing the color for the highest value.
#' @param mid A string informing the color for the mid value.
#' @param midpoint The mid value.
#' @param limit A vector containing c(lowest_value, highest_value).
#' @param name The matrix name.
#' 
#' @return The corresponding heat map.
matrix_vizualisation <- function(matrix_to_display, Mc, low, high, mid, midpoint, limit, name) {
  colnames(matrix_to_display) = colnames(Mc); rownames(matrix_to_display) = colnames(Mc)
  
  # Reorder the matrix from the biggest community to the smallest.
  communities = colnames(Mc[,order(-colSums(Mc))])
  matrix_to_display = round(matrix_to_display[communities, communities],3)
  
  # Create the chart.
  df_matrix_to_display = as.data.frame(as.table(as.matrix(matrix_to_display)))
  heat = ggplot(data = df_matrix_to_display, aes(x=Var1, y=Var2, fill=Freq)) + 
    geom_tile() +
    geom_text(aes(Var1, Var2, label = Freq), color = "black", size = 4) +
    scale_fill_gradient2(low = low, high = high, mid = mid, 
                         midpoint = midpoint, limit = limit,
                         name=name) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1)) +
    labs(title=name, x ="Community (from)", y = "Community (to)")
  
  print(heat)
  return (heat)
}

#' Create a cloud point chart to see the porosity of boundaries according to their size.
#' 
#' @param matrix_porosity The porosity matrix.
#' @param matrix_boundary_sizes The boundary sizes matrix.
#' 
#' @return The corresponding point cloud chart.
porosity_vizualisation2D <- function(matrix_porosity, matrix_boundary_sizes) {
  # Join the matrices.
  b_sizes = as.data.frame(as.table(matrix_boundary_sizes)); names(b_sizes)[names(b_sizes) == "Freq"] = "Size"
  p_scores = as.data.frame(as.table(as.matrix(matrix_porosity))); names(p_scores)[names(p_scores) == "Freq"] = "Porosity"
  df_scores = merge(b_sizes, p_scores, by=c("Var1", "Var2"), all.x=TRUE)
  names(df_scores)[names(df_scores) == "Var1"] = "Community"
  
  # Create the chart.
  pointcloud = ggplot(df_scores, aes(x=Size, y=Porosity, color=Community)) + geom_point() + labs(title="Size and porosity of boundaries", x ="Size (vertex count)", y = "Porosity (in %)") + stat_ellipse()
  
  print(pointcloud)
  return (pointcloud)
}