import matplotlib as mlab
mlab.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from benchmark import *

class BenchmarkPlot(object):
    def __init__(self, filename):
        self.filename = filename
        
    def plot_accuracy(self, lab, recall=0, precision=0):
        self.lab = lab
        self.recall = recall
        self.precision = precision

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

    def plot_editing_type(self):
        pass


    def plot_ag_enrichment(self):
        pass

        
    def plot_all_stats(self):
        pass
