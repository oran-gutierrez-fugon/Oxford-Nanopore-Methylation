# Load necessary libraries
library(ggplot2)
library(scales)
library(tools)

# New data frame with updated terms and p-values
data <- data.frame(
  Term = c("neuron projection (GO:0043005)", "axon (GO:0030424)",
           "asymmetric synapse (GO:0032279)", "postsynaptic density (GO:0014069)",
           "dendrite (GO:0030425)", "calcium channel complex (GO:0034704)",
           "glutamatergic synapse (GO:0098978)", "cytoskeleton of presynaptic active zone (GO:0048788)",
           "cilium (GO:0005929)", "cation channel complex (GO:0034703)"),
  P_Value = c(1.81E-17, 7.23E-16, 2.55E-10, 7.16E-10,
              1.28E-07, 1.23E-05, 2.12E-05, 0.000101,
              0.000202, 0.000241)

)

# Sort the data by P_Value in ascending order
data <- data[order(data$P_Value), ]

# Convert term names to title case
data$Term <- tools::toTitleCase(data$Term)

# Transform the p-values into a -log10 scale for better visualization
data$NegLog10_P_Value <- -log10(data$P_Value)

# Generate labels for the bars with the format: Term * p-value
data$labels <- paste(data$Term, "*", formatC(data$P_Value, format = "e", digits = 2), sep = " ")

# Create the bar plot with the new data and updated title
p <- ggplot(data, aes(x=NegLog10_P_Value, y=factor(Term, levels=rev(Term)))) +
  geom_bar(stat="identity", fill="lightgray") +
  geom_text(aes(x = 0.15, label = labels), hjust=0, size=4.2, color="black", position=position_stack(vjust=0.5)) +
  geom_segment(x = -0.02, y = 0.43, xend = max(data$NegLog10_P_Value) + 0.5, yend = 0.43, color = "darkgray", size = 0.5) +
  geom_segment(x = 0, y = 0.43, xend = 0, yend = length(unique(data$Term)) + 0.45, color = "darkgray", size = 0.5) + # Vertical line
  theme_minimal() +
  labs(x=expression("-log"[10]*"(p-value)"), y="", title="GO Cellular Component 2021 LUHMES Neurons") +
  theme(
    axis.text.x = element_text(color = "black", margin = margin(t = 2), size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_line(color = "darkgray", size = 0.5),
    axis.ticks.x.top = element_blank()
  ) +
  scale_x_continuous(
    breaks = seq(0, 16, by = 2),
    limits = c(0, max(data$NegLog10_P_Value) + 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off")


# Save the plot as a PDF to the specified path
pdf_path <- "GO_Cellular_Component_2021_LUHMES_Neurons.pdf"
ggsave(pdf_path, plot = p, width = 6, height = 6)

pdf_path
