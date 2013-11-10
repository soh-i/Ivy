library("ggplot2")
library("reshape2")

data <- read.table("ivy_bench.log", header=T, sep=",")

g <- ggplot(
  data,
  aes(x=Precision, y=Recall)) +
  geom_point(aes(colour=CSV), size=10, shape=19) +
  ylim(0,1) +
  xlim(0,1) +
  labs(
    title=paste("Benchmarking test for", "Human", "data set"),
    x="Precision",
    y="Recall"
  ) + 
  #theme_bw(base_size=20, base_family="Helvetica")# +
  #coord_fixed(ratio=1) +
  theme(axis.line.x=element_line(size=1),
        axis.line.y=element_line(size=1),
        axis.ticks.x=element_line(size=1),
        axis.ticks.y=element_line(size=1),
        axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18)
  #      legend.key=element_rect(colour="transparent", fill="transparent"),
  #      plot.background=element_rect(fill="transparent", colour="transparent"),
  #      legend.background=element_rect(fill="transparent", colour="transparent"),
  #      panel.background=element_rect(fill="transparent", colour="transparent")
  )
plot(g)

#ggsave(filename=(paste("benchmarking_human.png")), plot=g, height=10, width=10)
