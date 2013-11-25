human <- melt(read.table("human_bench.log", header=T, sep=","))
mouse <- melt(read.table("mouse_bench.log", header=T, sep=","))
fly <- melt(read.table("fly_bench.log", header=T, sep=","))

human_f_score <- subset(human, human$variable == "F.measure")
mouse_f_score <- subset(mouse, mouse$variable == "F.measure")
fly_f_score <- subset(fly, fly$variable == "F.measure")

bind_data <- rbind(human_f_score, mouse_f_score, fly_f_score)

order <- seq(length(bind_data$value), 1, by=-1)

g <- ggplot(
  bind_data,
  aes(y=sort(bind_data$value), x=order)) +
  geom_line(aes(colour=Species), size=2) +
  geom_point(colour="gray") + 
  labs(x="Rank order", y="F measure") + 
  ylim(0,1) + 
  theme(axis.line=element_line(size=1), # line weight
        axis.ticks=element_line(size=1, colour="black"), #x/y ticks size
        axis.text=element_text(size=20, colour="black"), # x/y line weight
        axis.title=element_text(size=20, colour="black"), #x/y lab font size
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        legend.key=element_rect(colour="transparent", fill="transparent"),
        plot.background=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="transparent"),
        panel.background=element_rect(fill="transparent", colour="transparent"))
plot(g)
ggsave(filename="rank_order.png", width=300, height=250, units="mm", plot=g, bg="transparent")
