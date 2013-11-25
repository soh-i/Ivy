library("ggplot2")
library("reshape2")

data <- read.table("ivy_bench.csv", header=T, sep=",")

g <- ggplot(
  data,
  aes(x=Precision, y=Recall)) +
  geom_point(aes(colour=CSV), size=9, shape=19) +
  stat_smooth(method=lm, se=TRUE, level=0.95) + 
  ylim(0,1) +
  xlim(0,1) +
  labs(
    #title=paste("Benchmarking test for", "Human", "data set"),
    x="Precision",
    y="Recall"
  ) + 
  #theme_bw(base_size=20, base_family="Helvetica") +
  coord_fixed(ratio=1) +
  theme(axis.line=element_line(size=1), # line weight
        axis.ticks=element_line(size=1, colour="black"), #x/y ticks size
        axis.text=element_text(size=18, colour="black"), # x/y line weight
        axis.title=element_text(size=18, colour="black"), #x/y lab font size
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        legend.title=element_text(size=18),
        legend.text=element_text(size=18),
        legend.key=element_rect(colour="transparent", fill="transparent"),
        plot.background=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="transparent"),
        panel.background=element_rect(fill="transparent", colour="transparent")
                  )
plot(g)

ggsave(filename="test.png", width=400, height=300, units="mm", plot=g, bg="transparent")

