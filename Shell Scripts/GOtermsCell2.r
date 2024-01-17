# Adjusted segment for creating the bar plot
p <- ggplot(data, aes(x=NegLog10_P_Value, y=Component)) +
  geom_bar(stat="identity", fill="lightgray") +
  geom_text(aes(x = 0.2, label = labels), hjust=0, size=3.5, color="black", position=position_stack(vjust=0.5)) +
  geom_segment(x = 0, y = 0.45, xend = 0, yend = length(unique(data$Component)) + 0.45, color = "darkgray") + # Vertical line
  geom_segment(x = 0, y = 0.45, xend = max(data$NegLog10_P_Value), yend = 0.45, color = "darkgray") + # Horizontal line
  theme_minimal() +
  labs(x=expression("-log"[10]*"(p-value)"), y="", title="GO Cellular Component 2021") +
  theme(
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_blank(),
    plot.title = element_text(hjust = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.ticks.x = element_line(color = "darkgray"),
    axis.ticks.x.top = element_blank()
  ) +
  scale_x_continuous(
    breaks = seq(0, ceiling(max(data$NegLog10_P_Value)), by = 2),
    expand = c(0, 0), # Remove space between bars and axis
    sec.axis = dup_axis(labels = NULL, name = NULL)
  ) +
  coord_cartesian(clip = "off") # Adjust plot area to avoid clipping of lines

# Print the plot
print(p)