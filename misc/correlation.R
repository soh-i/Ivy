library("ggplot2")
library("reshape2")

#data <- melt(read.table("fly_bench.log", header=T, sep=","))
data <- melt(read.table("human_bench.log", header=T, sep=","))
r <- subset(data, data$variable == "Recall")
p <- subset(data, data$variable == "Precision")
c <- subset(data, data$variable == "PredCount")
d <- rbind(p, r)
PredCount <- c$value
d <- cbind(d, PredCount)

names(d)[names(d) == "CSV"] <- "Study"
names(d)[names(d) == "variable"] <- "Metrics"

g <- ggplot(d, aes(x=value, y=PredCount)) +
  geom_point(aes(colour=Study, shape=Metrics), size=10)  + 
  stat_smooth(method="lm", se=F, fullrange=F, aes(fill=Metrics)) + 
  labs(x="Metrics score", y="Candidate data count") +
  theme(axis.line=element_line(size=1), # line weight
        axis.ticks=element_line(size=1, colour="black"), #x/y ticks size
        axis.text=element_text(size=24, colour="black"), # x/y line weight
        axis.title=element_text(size=24, colour="black"), #x/y lab font size
        panel.grid.major=element_line(colour=NA),
        panel.grid.minor=element_line(colour=NA),
        legend.title=element_text(size=24),
        legend.text=element_text(size=24),
        legend.key=element_rect(colour="transparent", fill="transparent"),
        plot.background=element_rect(fill="transparent", colour="transparent"),
        legend.background=element_rect(fill="transparent", colour="transparent"),
        panel.background=element_rect(fill="transparent", colour="transparent"))

#ggsave(filename="correlations.png", width=440, height=300, units="mm", plot=g, bg="transparent")

#in human data
ggsave(filename="correlations.png", width=620, height=360, units="mm", plot=g, bg="transparent")

plot(g)

