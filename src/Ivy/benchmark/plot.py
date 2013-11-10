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
        axes.patch.set_facecolor('0.85')
        axes.set_axisbelow(True)

        plt.grid(True, which='minor', color="0.92", linestyle="-", linewidth=0.7)
        plt.grid(True, which='major', color="w", linestyle="-", linewidth=1.2)
        
        plt.xlabel("Precision")
        plt.ylabel("Recall")
        plt.ylim(0.0,1.0)
        plt.xlim(0.0,1.0)
        
        # remove axis border line
        for child in axes.get_children():
            # TODO: do not use mlab
            import matplotlib
            if isinstance(child, matplotlib.spines.Spine):
                child.set_alpha(0)
                
        # background
        for line in axes.get_xticklines() + axes.get_yticklines():
            line.set_markersize(5)
            line.set_color("gray")
            line.set_markeredgewidth(1.4)
        
        plt.plot(self.precision, self.recall,
                 marker='o', color="red",
                 linestyle=".", markersize=10,
                 markeredgecolor="red", label=self.lab)
        plt.title("Benchmarking test for detection accuracy in " + self.sp)
        
        plt.legend()
        if axes.legend_ is not None:
            lg = axes.legend_
            lg.get_frame().set_linewidth(0)
            lg.get_frame().set_alpha(0.5)

        fig.savefig(self.filename + '.pdf')

    def plot_editing_type(self):
        pass

    def plot_ag_enrichment(self):
        pass
        
    def plot_all_stats(self):
        pass

if __name__ == '__main__':
    bplot = BenchmarkPlot("test", "human")
    bplot.plot_accuracy(lab=["lab1", "lab2"], recall=[0.88, 0.21], precision=[0.31, 0.78])
