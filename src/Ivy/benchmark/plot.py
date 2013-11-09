import matplotlib as mlab
mlab.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

__program__ = 'ivy_benchmark'
__author__ = 'Soh Ishiguro <yukke@g-language.org>'
__license__ = ''
__status__ = 'development'


class BenchmarkPlot(object):
    '''
    BenchmarkPlot class is to vizualize some stats of benchmarking results.
    >>> bplt = BenchmarkPlot("out_file")
    '''
    
    def __init__(self, filename, sp):
        self.filename = filename
        self.sp = sp
        
    def plot_accuracy(self, lab=None, recall=0, precision=0):
        '''
        >>> bplt.plot_accuracy(lab="sample01", recall=0.8, precision=0.78)
        Generate out_file.pdf
        '''
        
        self.lab = lab
        self.recall = recall
        self.precision = precision
        fig = plt.figure()
        
        axes = plt.subplot(111)
        axes.spines['right'].set_color('none')
        axes.spines['top'].set_color('none')
        axes.xaxis.set_ticks_position('bottom')
        axes.yaxis.set_ticks_position('left')
        
        axes.xaxis.grid(False)
        axes.yaxis.grid(False)
        
        plt.plot(self.precision, self.recall,
                 marker='o', color="black",
                 linestyle=".", markersize=10,
        markeredgecolor="black", label=self.lab)
        plt.title("Benchmarking test for detection accuracy in " + self.sp)
        
        plt.xlabel("Precision")
        plt.ylabel("Recall")
        plt.ylim(0.0,1.0)
        plt.xlim(0.0,1.0)
        plt.grid(color="gray")
        plt.legend()
        fig.savefig(self.filename + '.pdf')

    def plot_editing_type(self):
        pass

    def plot_ag_enrichment(self):
        pass
        
    def plot_all_stats(self):
        pass

