
## center ggplot title
.center_title = function() {
  ggplot2::theme(plot.title = element_text(hjust = 0.5))
}

.rotate_x_labels = function() {
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
}
