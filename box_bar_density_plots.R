
# --------- Expression labels ---------
display.labels <- c(
  "a" = expression(a),
  "rx" = expression(r[X]),
  "theta_X" = expression(theta[X]),
  "sigma_Y" = expression(sigma[Y])
)

# --------- Bar plot function ---------
make_barplot <- function(data, colnames_vec, ylab_expr, ymax) {
  if (is.vector(data)) {
    data <- as.data.frame(t(data))
    colnames(data) <- colnames_vec
  } else if (is.matrix(data) && is.null(colnames(data))) {
    colnames(data) <- colnames_vec
  }
  
  df <- tidyr::pivot_longer(as.data.frame(data), cols = everything(), names_to = "label", values_to = "value")
  df$label <- factor(df$label, levels = colnames_vec)
  
  ggplot(df, aes(x = label, y = value)) +
    geom_bar(stat = "summary", fun = mean, fill = "gray", color = "black", width = 0.8) +
    scale_x_discrete(labels = display.labels) +
    coord_cartesian(ylim = c(0, ymax)) +
    ylab(ylab_expr) +
    xlab(NULL) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# --------- Boxplot function ---------
make_boxplot <- function(data, colnames_vec, ylab_expr, ylim_vals = NULL) {
  data <- matrix(data, nrow = nrow(as.matrix(data)))
  data <- data[, seq_along(colnames_vec)]
  colnames(data) <- colnames_vec
  
  df <- pivot_longer(as.data.frame(data), cols = everything(), names_to = "label", values_to = "value")
  df$label <- factor(df$label, levels = colnames_vec)
  
  p <- ggplot(df, aes(x = label, y = value)) +
    geom_boxplot(fill = "lightgray", color = "black") +
    scale_x_discrete(labels = display.labels) +
    ylab(ylab_expr) +
    xlab(NULL) +
    theme_minimal(base_size = 16) +
    theme(
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16),
      plot.margin = margin(5, 5, 5, 5)
    )
  if (!is.null(ylim_vals)) {
    p <- p + coord_cartesian(ylim = ylim_vals)
  }
  return(p)
}

# --------- Density plot function ---------
make_densityplot <- function(data, colnames_vec, ylab_expr = NULL, ymax = NULL) {
  error_vector <- as.vector(data)
  
  ggplot(data.frame(error = error_vector), aes(x = error)) +
    geom_density(color = "red", linewidth = 0.5) +
    labs(x = "Error", y = NULL, title = ylab_expr) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title.y = element_blank(),
      panel.grid = element_blank()
    )
}