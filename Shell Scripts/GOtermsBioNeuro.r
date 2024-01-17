# Load necessary libraries
library(ggplot2)
library(scales)
library(tools)

# New data frame with updated terms and p-values
data <- data.frame(
  Term = c("nervous system development (GO:0007399)", "axonogenesis (GO:0007409)",
           "synapse organization (GO:0050808)", "axon guidance (GO:0007411)",
           "modulation of chemical synaptic transmission (GO:0050804)", 
           "chemical synaptic transmission (GO:0007268)", "positive regulation of axonogenesis (GO:0050772)",
           "neuron projection development (GO:0031175)", "synapse assembly (GO:0007416)",
           "neuron projection morphogenesis (GO:0048812)"),
  P_Value = c(1.62E-17, 2.05E-15, 3.08E-15, 1.24E-13, 1.00E-11, 5.09E-10, 8.24E-10, 2.92E-09, 3.19E-09, 4.57E-09)
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
  labs(x=expression("-log"[10]*"(p-value)"), y="", title="GO Biological Process 2021 Upregulated in LUHMES Neurons") +
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
pdf_path <- "GO_Biological_Process_2021_Upregulated_in_LUHMES_Neurons.pdf"
ggsave(pdf_path, plot = p, width = 6, height = 6)

pdf_path
