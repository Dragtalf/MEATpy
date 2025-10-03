# Patricia Tran (METABOLIC)
# Updated by BMF 2025-04-25 (MEATpy)
# Generate a network plot with reactions as nodes and MAGs as edges, colored by taxonomic group.

userprefs <- commandArgs(trailingOnly = TRUE)
R_input_table <- userprefs[1]
plots.folder.path <- userprefs[2]

if (length(userprefs) > 2){
  mirror.location <- userprefs[3]
} else {
  mirror.location <- "https://cran.mtu.edu"
}

network.plots.folder <- file.path(plots.folder.path, "network_plot")

# Clean folder if it already exists
if (dir.exists(network.plots.folder)) {
  unlink(network.plots.folder, recursive = TRUE)
}
dir.create(network.plots.folder, recursive = TRUE)

# Load libraries silently
suppressPackageStartupMessages({
  library(ggraph)
  library(igraph)
  library(tidyverse)
  library(tidygraph)
})

# Read input
table <- read.csv(R_input_table, header = TRUE, sep = "\t")

# Build graph
my_graph <- table[, c(2, 3, 4, 5)] %>% graph_from_data_frame()
deg <- degree(my_graph, mode = "all")

# Whole-community plot
community.plot <- table[, c(2, 3, 4, 5)] %>% 
  graph_from_data_frame() %>% 
  ggraph(layout = "linear", circular = TRUE) +
  geom_edge_arc(alpha = 0.25, aes(width = Coverage.value.average., color = as.factor(Taxonomic.Group))) +
  geom_node_point(aes(size = 0.02 * deg), color = "black", alpha = 0.75) +
  geom_node_text(aes(label = name), color = "black", repel = TRUE) +
  scale_edge_width_continuous(trans = "identity") +
  scale_edge_color_discrete() +
  labs(title = "Metabolic connections within dataset", subtitle = "No scaling") +
  theme_graph()

# Save community plot
plot.name <- file.path(network.plots.folder, "CommunityPlot.PDF")
if (file.exists(plot.name)) file.remove(plot.name)
cairo_pdf(filename = plot.name, width = 11, height = 8.5, onefile = TRUE)
print(community.plot)
dev.off()

# Individual taxonomic group plots
taxo.groups <- unique(table$Taxonomic.Group)

for (taxo in taxo.groups) {
  message("Making a plot for: ", taxo)
  
  individual.table <- subset(table, Taxonomic.Group == taxo)
  
  ind.plot <- individual.table[, c(2, 3, 4, 5)] %>%
    graph_from_data_frame() %>%
    ggraph(layout = "fr") +
    geom_edge_link(alpha = 0.5, aes(width = Coverage.value.average.)) +
    geom_node_point(color = "black", size = 2) +
    geom_node_text(aes(label = name), color = "black", repel = TRUE) +
    theme_graph() +
    labs(
      title = paste0("Functional connections within ", taxo),
      subtitle = "No scaling"
    )
  
  plot.name2 <- file.path(network.plots.folder, paste0(taxo, ".Individual.Taxonomic.Groups.Functional.Network.PDF"))
  if (file.exists(plot.name2)) file.remove(plot.name2)
  
  cairo_pdf(filename = plot.name2, width = 11, height = 8.5, onefile = TRUE)
  print(ind.plot)
  dev.off()
}
