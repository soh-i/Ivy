import matplotlib as mlab
mlab.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import (
    MultipleLocator,
    FormatStrFormatter,
    )

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
        
        plt.xlabel("Precision")
        plt.ylabel("Recall")
        plt.ylim(0.0,1.0)
        plt.xlim(0.0,1.0)
        
        axes.patch.set_facecolor('0.89')
        axes.set_axisbelow(True)
        
        axes.xaxis.set_major_locator(MultipleLocator(0.2))
        axes.xaxis.set_minor_locator(MultipleLocator(0.1))
        axes.yaxis.set_major_locator(MultipleLocator(0.2))
        axes.yaxis.set_minor_locator(MultipleLocator(0.1))
        axes.xaxis.grid(True, "minor")
        axes.yaxis.grid(True, "minor")
        axes.grid(b=True, which='major', color="w", linestyle="-", linewidth=1.4)
        axes.grid(b=True, which='minor', color="0.97", linestyle="-", linewidth=0.5)
                                     
        for line in axes.xaxis.get_ticklines(minor=True) + axes.yaxis.get_ticklines(minor=True):
            line.set_markersize(0)

        # remove axis border line
        for child in axes.get_children():
            # TODO: do not use import statement
            import matplotlib
            if isinstance(child, matplotlib.spines.Spine):
                child.set_alpha(0)
                
        # background
        for line in axes.get_xticklines() + axes.get_yticklines():
            line.set_markersize(5)
            line.set_color("gray")
            line.set_markeredgewidth(1.4)

        # generate color map
        colors = []
        for i in range(len(self.precision)):
            colors.append(cm.PiYG(float(i)/len(self.precision),1))
            
        for _ in range(len(self.precision)):
            plt.plot(self.precision[_], self.recall[_],
                     color=colors[_],
                     marker='o',
                     linestyle=".", markersize=10,
                     markeredgecolor=colors[_],
                     label=self.lab[_])
        
        plt.title("Benchmarking test for detection accuracy in " + self.sp)
        
        plt.legend(fontsize=10)
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

    
