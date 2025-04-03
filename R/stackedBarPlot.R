library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)

data <- read_excel('rutils/microvessel_celltype.xlsx')

data_long <- pivot_longer(data, cols = day0:day35, 
                          names_to = "Day", 
                          values_to = "Percent") %>%
  mutate(Condition = ifelse(CellType == 'Proliferating', "Proliferating", "Not-Proliferating")) # %>%
  #group_by(Day, Condition) %>%
  #summarize(TotalPercent = sum(Percent, na.rm = TRUE))

bPlot <- ggplot(data_long, 
                aes(fill=Condition, # rearrange
                    y=Day,
                    x=Percent)) + 
  scale_fill_manual(values = c("Proliferating" = "red", "Not-Proliferating" = "grey")) +
  ylab("Day") +
  xlab("CellType Ratios") +
  geom_bar(position="fill", stat="identity", color = 'black') + 
  labs(fill = "Cluster") +
  theme_classic() + 
  theme(axis.title.x = element_text(size = 8, face = "bold", colour = 'black'),
        axis.text.x = element_text(size = 8, face = "bold", colour = 'black'),
        axis.title.y = element_text(size=8, face="bold"),
        axis.text.y = element_text(size = 8, face = "bold", colour = 'black'),
        legend.title = element_blank(),
        legend.text = element_text(size=8),
        legend.position = 'none')

ggsave(bPlot, 
       filename = "rutils/angela_plot_lines.svg", 
       width = 8, 
       height = 4)
