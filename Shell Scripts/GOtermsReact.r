# Load necessary libraries
library(ggplot2)
library(scales)
library(tools)

# Data frame containing the cellular components and their corresponding p-values
data <- data.frame(
  Term = c("Neuronal System (R-HSA-112316)", "Transmission Across Chemical Synapses (R-HSA-112315)",
           "Dopamine Neurotransmitter Release Cycle (R-HSA-212676)", "Neurexins And Neuroligins (R-HSA-6794361)",
           "Protein-protein Interactions At Synapses (R-HSA-6794362)", "Serotonin Neurotransmitter Release Cycle (R-HSA-181429)",
           "Acetylcholine Neurotransmitter Release Cycle (R-HSA-264642)", "Glutamate Neurotransmitter Release Cycle (R-HSA-210500)",
           "Phase 0 - Rapid Depolarisation (R-HSA-5576892)", "CRMPs In Sema3A Signaling (R-HSA-399956)"),
  P_Value = c(1.77e-14, 2.64e-09, 2.79e-09, 1.04e-08, 3.67e-08, 6.30e-08, 2.93e-06, 3.48e-06, 5.18e-06, 9.05e-06)
)

# Sort the data by P_Value in ascending order
data <- data[order(data$P_Value), ]

# Transform the p-values into a -log10 scale for better visualization
data$NegLog10_P_Value <- -log10(data$P_Value)

# Generate labels for the bars with the format: Term * p-value
data$labels <- paste(data$Term, "*", formatC(data$P_Value, format = "e", digits = 2), sep = " ")

# Create the bar plot
p <- ggplot(data, aes(x=NegLog10_P_Value, y=factor(Term, levels=rev(Term)))) +
  geom_bar(stat="identity", fill="lightgray") +
  geom_text(aes(x = 0.15, label = labels), hjust=0, size=4.2, color="black", position=position_stack(vjust=0.5)) +
  geom_segment(x = -0.02, y = 0.43, xend = max(data$NegLog10_P_Value) + 0.5, yend = 0.43, color = "darkgray", size = 0.5) + # Horizontal line
  geom_segment(x = 0, y = 0.43, xend = 0, yend= length(unique(data$Term)) + 0.45, color = "darkgray", size = 0.5) + # Vertical line
  theme_minimal() +
  labs(x=expression("-log"[10]*"(p-value)"), y="", title="Reactome 2022 LUHMES Neurons") +
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
    breaks = seq(0, 14, by = 2),
    limits = c(0, max(data$NegLog10_P_Value) + 1),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off")
               

# Save the plot as a PDF to the specified path
pdf_path <- "Reactome_2022_LUHMES_Neurons.pdf"
ggsave(pdf_path, plot = p, width = 6, height = 6)

pdf_path

