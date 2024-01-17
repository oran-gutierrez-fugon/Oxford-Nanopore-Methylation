# Load necessary libraries
library(ggplot2)
library(scales)
library(tools)

# New data frame with updated terms and p-values
data <- data.frame(
  Term = c("ribosome biogenesis (GO:0042254)", "gene expression (GO:0010467)",
           "translation (GO:0006412)", "rRNA processing (GO:0006364)",
           "cellular macromolecule biosynthetic process (GO:0034645)",
           "ncRNA processing (GO:0034470)", "mRNA processing (GO:0006397)",
           "mRNA splicing, via spliceosome (GO:0000398)",
           "RNA Splicing via Bulged Ade Transesterification (GO:0000377)",
           "rRNA metabolic process (GO:0016072)"),
  P_Value = c(2.37E-76, 1.19E-72, 2.47E-72, 3.93E-72,
              1.75E-71, 3.86E-65, 7.55E-65, 9.47E-65,
              3.49E-62, 1.13E-60)
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
  labs(x=expression("-log"[10]*"(p-value)"), y="", title = "GO Biological Process 2021 Downregulated in LUHMES Neurons") +
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
    breaks = seq(0, 72, by = 8),
    limits = c(0, max(data$NegLog10_P_Value) + 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off")


# Save the plot as a PDF to the specified path
pdf_path <- "GO_Biological_Process_2021_Downregulated_in_LUHMES_Neurons.pdf"
ggsave(pdf_path, plot = p, width = 6, height = 6)

pdf_path
