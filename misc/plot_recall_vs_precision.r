library("ggplot2")
library("reshape2")

#data <- read.table("/Volumes/TTCK.cs0/benchmark/results/bench_DARNED_mouse.log", header=T)
sp <- ""

g <- ggplot(
  lympho_peng,
  aes(x=Precision, y=Recall)) +
  #geom_point(aes(colour=Label), size=8, shape=19) +
  geom_point(aes(colour=AnswerSet), size=4, shape=19) +
  ylim(0,1) +
  xlim(0,1) +
  labs(#title=paste("Benchmarking test for", sp, "data set"),# sp, "data sets"),
    x="Precision",
    y="Recall"
  ) + #scale_colour_brewer(palette="Set2") + 
  theme_classic(base_size=20, base_family="Helvetica") +
  coord_fixed(ratio=1) +
  theme(axis.line.x=element_line(size=1),
        axis.line.y=element_line(size=1),
        axis.ticks.x=element_line(size=1),
        axis.ticks.y=element_line(size=1),
        axis.text.y=element_text(size=18),
        axis.text.x=element_text(size=18),
        legend.key=element_rect(colour="transparent", fill="transparent"),
        plot.background=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="transparent"),
        panel.background=element_rect(fill="transparent", colour="transparent")
  )

plot(g)

#ggsave(filename=(paste("benchmarking_human.png")), plot=g, height=10, width=10)

