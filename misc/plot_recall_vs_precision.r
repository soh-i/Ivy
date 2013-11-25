library("ggplot2")
library("reshape2")

plot_benchmarking <- function(data, filename) {
  label = "Study"
  names(data)[names(data) == "CSV"] <- label
  
  g <- ggplot(
    data,
    aes(x=Precision, y=Recall)) +
    geom_point(aes(colour=Study), size=9, shape=19) +
    #stat_smooth(method=lm, se=FALSE) +
    ylim(0,1) +
    xlim(0,1) +
    labs(
      #title=paste("Benchmarking test for", "Human", "data set"),
      x="Precision",
      y="Recall"
    ) + 
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
  ggsave(filename=paste(filename, ".png", sep=""), width=400, height=300, units="mm", plot=g, bg="transparent")
}

plot_all_sp <- function() {
  human_data <- read.table("ivy_bench.csv", header=T, sep=",")
  plot_benchmarking(human_data, "bench_human")
  
  mouse_data <- read.table("../mouse_bench.log", header=T, sep=",")
  plot_benchmarking(mouse_data, "bench_mouse")

  fly_data <- read.table("../fly_bench.log", header=T, sep=",")
  plot_benchmarking(fly_data, "bench_fly")
}

plot_all_sp()


