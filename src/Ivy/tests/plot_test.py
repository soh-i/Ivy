from Ivy.benchmark.benchmark import *

plot = BenchmarkPlot("test_data", recall=[0.78, 0.92], precision=[0.12,0.48])
plot.p_r_plot()
