import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plot
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties
import scipy
import scipy.cluster.hierarchy as sch
from operator import *
pd.core.config.option_context('mode.use_inf_as_null',True)
try:
    import seaborn as sns  
    rcdefsns = plot.rcParams.copy()
except:
    pd.options.display.mpl_style = 'default'
    plot.ioff()


def mhist(values,*args):
    return np.histogram(values,args[0],normed=True)[0]
def binning(maxlength,size):    
    bins=np.arange(5,maxlength,size)
    return bins

def table(filein,metric,counts,maxlength,size):
    table = pd.io.parsers.read_table(filein,index_col=0)
    metric='A_length_'+str(metric)
    table[metric]=table[metric].str.split(',').apply(np.asarray,dtype=int)    
    index=table['No.Aligned_Reads'].where(table['No.Aligned_Reads']>=counts).dropna().index
    table=table.ix[index]
    return table

def common_index(tables):
    index=tables[0][1].index
    for i in range(1,len(tables)):
        index=index.intersection(tables[i][1].index)
    return index

def ksstat(control,experimental):
    return scipy.stats.ks_2samp(control,experimental)

def ks(tables,metric,control,minks,pvalue):
    metric='A_length_'+str(metric)
    index=common_index(tables)
    ks_data=pd.DataFrame(index=index)
    for i in tables:
        ks_data[i[0]]=i[1].ix[index][metric]

    for column in ks_data.columns:
        ks='KS-Stat'+str(control)+str(column)
        pval='Pvalue'+str(control)+str(column)

        ks_data[ks],ks_data[pval]=zip(*ks_data.apply(lambda x:ksstat(x[control],x[column]),axis=1))

        ks_data[ks]=ks_data[ks].where(ks_data[ks]>=minks,0.0)
        ks_data[ks]=ks_data[ks].where(ks_data[pval]<=pvalue,0.0)

    return ks_data        

def heatmap(tables,labels,maxlength,size):
    fig, axes = plot.subplots(nrows=1, ncols=len(labels)+1)

    for i in range(len(labels)):
        datap=pd.DataFrame(index=tables.index)
        datap[labels[i]]=tables[labels[i]]
        datap['HIST']=datap[labels[i]].apply(mhist,args=(binning(maxlength,size),))
        dataplot=np.row_stack(datap['HIST'].values)
        heatmap = axes[i].matshow(dataplot, aspect='auto', cmap=plot.cm.Greens,vmin=0,vmax=0.05)
        axes[i].set_xticks([])
        axes[i].set_yticks([])
        axes[i].set_title(str(labels[i]))
        
    ksplot=np.row_stack(tables.filter(regex='KS').values)    
    heatmap = axes[-1].matshow(ksplot, aspect='auto', cmap=plot.cm.Reds,vmin=0,vmax=1.)
    axes[-1].set_xticks([])
    axes[-1].set_yticks([])
    axes[-1].set_title(str('KS-STAT'))
    return fig
    

def control_sample(labels,control):
    if control == '':
        return labels[0]
    else:
        if control in labels:
            return control
        else:
            print 'Control Label is not present in Labels revert to sample %s' %(options.label[0])
            return label[0]    
        
def sort_on_ks(tables,control,sort_sample):
    column='KS-Stat'+str(control)+str(sort_sample)
    tables.sort([column],ascending=True,inplace=True)
    return tables

def run(parser,options):
    data=[]

    for i in range(len(options.tables[0])):
        data.append((options.labels[0][i],
                       table(options.tables[0][i],
                             options.metric,
                             options.counts,
                             options.max_length,
                             options.binsize)))
    control=control_sample(options.labels[0],options.control)
    ksdata=ks(data,options.metric,control,options.minks,options.pvalue)
    sortdata=sort_on_ks(ksdata,control,options.sort)

    fig_pdf=PdfPages('%s.pdf' %('KS_heatmap'))
    fig=heatmap(sortdata,
                options.labels[0],
                options.max_length,
                options.binsize)
    fig_pdf.savefig()
    plot.close(fig)
    fig_pdf.close()        
