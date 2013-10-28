import matplotlib as mlab
mlab.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from benchmark import *

class BenchmarkPlot(object):
    def __init__(self, lab, recall=None, precision=None):
        self.recall = recall
        self.precision = precision
        self.lab = lab
        
    def p_r_plot(self):
        fig = plt.figure()
        
        plt.plot(self.precision, self.recall,
                 marker='o', color="green",
                 linestyle=".", markersize=10,
                 markeredgecolor="green", label=self.lab)
        plt.title("Benchmarking test")
        plt.xlabel("Precision")
        plt.ylabel("Recall")
        plt.ylim(0.0,1.0)
        plt.xlim(0.0,1.0)
        plt.legend()
        fig.savefig("test.pdf")

if __name__ == '__main__':
    plot = BenchmarkPlot("test_data", recall=[0.78, 0.92], precision=[0.12,0.48])
    plot.p_r_plot()
