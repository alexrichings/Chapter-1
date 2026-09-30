
# Chapter 1 

#------------------------------------------------------------------------------
  
# R: plot_helpers ---- 

# Purpose: 
# Help with adding standard error ratio (ratio of standard errors between two odds ratios) to a column in a forest plot 

#------------------------------------------------------------------------------

# function 
# add SE column to forest plot 
  
add_se_column <- function(plot, se_labels,
                          right_col = 7, strip_row = 8,
                          top_row = 10, bot_row = 10,
                          col_width = 2.8) {
  gt     <- ggplotGrob(plot)
  n_rows <- length(se_labels)
  
  header_grob <- grobTree(
    rectGrob(gp = gpar(fill = "grey85", col = "grey35", lwd = 0.8 * .pt)),
    textGrob("SE ratio", gp = gpar(fontface = "bold", fontsize = 18))
  )
  
  row_grobs <- lapply(se_labels, function(lab) {
    grobTree(
      rectGrob(gp = gpar(fill = "white", col = "grey35", lwd = 0.8 * .pt)),
      textGrob(lab, gp = gpar(fontsize = 18), x = 0.5, y = 0.5)
    )
  })
  
  se_column_grob <- frameGrob(
    layout = grid.layout(nrow = n_rows, ncol = 1,
                         heights = unit(rep(1, n_rows), "null"))
  )
  for (i in seq_along(row_grobs)) {
    se_column_grob <- placeGrob(se_column_grob, row_grobs[[i]], row = i, col = 1)
  }
  
  gt      <- gtable_add_cols(gt, unit(col_width, "cm"), pos = right_col)
  new_col <- right_col + 1
  gt <- gtable_add_grob(gt, header_grob, t = strip_row, b = strip_row,
                        l = new_col, r = new_col, name = "se_header")
  gt <- gtable_add_grob(gt, se_column_grob, t = top_row, b = bot_row,
                        l = new_col, r = new_col, name = "se_cells")
  gt
}



