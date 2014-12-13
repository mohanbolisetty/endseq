'''endseq_main.py'''
import sys
import argparse

import endseq.version

def run_subtool(parser, args):
    if args.command == 'analysis':
        import analysis as submodule
    elif args.command == 'plots':
        if len(args.tables[0]) == 0:
            parser.print_help()
            parser.error('Requires atleast 1 data table.')
        if len(args.labels[0]) == 0:
            parser.print_help()
            parser.error('Requires atleast 1 label.')
        if len(args.tables[0]) != len(args.labels[0]) :
            parser.print_help()
            parser.error('No. of tables and labels are not the same')

        if args.plottype=='SingleGene':
            if args.geneid == '':
                parser.print_help()
                parser.error('Single gene plot requires GENE ID')
            else:
                args.counts=1
        
        import plots as submodule
    elif args.command == 'stats':
        import stats as submodule

    submodule.run(parser, args)

class ArgumentParserWithDefaults(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(ArgumentParserWithDefaults, self).__init__(*args, **kwargs)
        self.add_argument("-q", "--quiet", help="Do not output warnings to stderr",
                        action="store_true",
                        dest="quiet")

def csv(values):
    return values.split(',')

def main():

    parser = argparse.ArgumentParser(prog='endseq', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-v", "--version", help="Installed endseq version",
                    action="version",
                    version="%(prog)s " + str(endseq.version.__version__))

    subparsers = parser.add_subparsers(title='[sub-commands]', dest='command',parser_class=ArgumentParserWithDefaults)

    ###########
    # ANALYSIS
    ###########

    parser_analysis=subparsers.add_parser('analysis',
                                          help='endseq Alignments and A-Tail counting')
    parser_analysis.add_argument('-1', '--read1',
                                 dest='read1',
                                 required=True,
                                 help='Read 1 Fastq file',
                                 default=None)
    parser_analysis.add_argument('-2', '--read2', dest='read2',
                                 required=True, help='Read 2 Fastq file',
                                 default=None)

    parser_analysis.add_argument('-i', '--index', dest='index',
                                 required=True, help='path to bowtie index',
                                 default=None)
    parser_analysis.add_argument('-G', '--GFF', dest='GFF',
                                 required=True, help='Gene Features File',
                                 default=None)
    parser_analysis.add_argument('-o', '--output', '--output-folder',
                                 dest='outputfolder',required=True,
                                 help='output folder',
                                 default=None)    

            #######################
            # ADDITIONAL PARAMETERS
            #######################

    parser_analysis.add_argument('-a', '--adapter_sequence', dest='adapter',
                                 help='3`end adapter sequence',
                                 default='GATGGTGCCTACAGTTTTT')
    parser_analysis.add_argument('-s', '--align', dest='r1len',
                                 help='Length of read 1 for alignment',
                                 default=24,type=int)

    parser_analysis.add_argument('-e', '--error', dest='error',
                      help='Error rate for A counting',
                      default=10,type=int)
    
    parser_analysis.add_argument('-t', '--threshold', dest='threshold',
                                 help='Threshold for A percent',
                                 default=4,type=int)
    parser_analysis.add_argument('-w', '--window', dest='window',
                                 help='Window Size for A percent',
                                 default=10,type=int)

    parser_analysis.add_argument('-S', '--string', dest='string',
                                 help='String Length for A tail measurement',
                                 default=5,type=int)
    parser_analysis.add_argument('-d', '--distance', dest='distance',
                                 help='Max distance between strings',
                                 default=20,type=int)

    parser_analysis.set_defaults(func=run_subtool)
    

    ###########
    # PLOTS
    ###########

    parser_plots=subparsers.add_parser('plots',
                                          help='Plots for endseq analysis')
    parser_plots.add_argument('-t', '--tables', dest='tables',
                              help='Comma seperated list of analysis_tables',
                              required=True, nargs='*',type=csv,
                              default=None)
    parser_plots.add_argument('-l', '--label', dest='labels',
                              help='Comma seperated list of lables ',
                              required=True, nargs='*',type=csv,
                              default=None)

    parser_plots.add_argument('-p', '--plottype', dest='plottype',
                              help='Cumulative plots ', default='CDF',
                              choices=['KDE',
                                       'CDF',
                                       'HIST',
                                       'Heatmap',
                                       'Genewise',
                                       'ScatterMatrix',
                                       'SingleGene'])
    parser_plots.add_argument('-m', '--metric', dest='metric',
                              help='A length metric ',
                              choices=['counts','window','strings'],
                              default='strings')
    parser_plots.add_argument('-a', '--max_length',dest='max_length',
                              help='Max length of A tail',
                              default=200,type=int)
    parser_plots.add_argument('-s','--binsize',dest='binsize',
                              help='Size of Bins',default=1,type=int)
    parser_plots.add_argument('-r','--readsthreshold',dest='counts',
                              help='Min read counts',default=30,type=int)
    parser_plots.add_argument('-i','--geneid',
                              dest='geneid',help='Gene ID',
                              default=None)
    parser_plots.set_defaults(func=run_subtool)


    ###########
    # STATS
    ###########

    parser_stats=subparsers.add_parser('stats',
                                          help='Plots for endseq analysis')
    parser_stats.add_argument('-t', '--tables', dest='tables',
                              help='Comma seperated list of analysis_tables',
                              required=True, nargs='*',type=csv,
                              default=None)
    parser_stats.add_argument('-l', '--label', dest='labels',
                              help='Comma seperated list of lables ',
                              required=True, nargs='*',type=csv,
                              default=None)

    parser_stats.add_argument('-c','--controlsample',
                              dest='control',help='Control sample')
    parser_stats.add_argument('-m', '--metric', dest='metric',
                              help='A length metric ',
                              choices=['counts','window','strings'],
                              default='strings')
    parser_stats.add_argument('-a', '--max_length',dest='max_length',
                              help='Max length of A tail',
                              default=200,type=int)
    parser_stats.add_argument('-s','--binsize',dest='binsize',
                              help='Size of Bins',default=5,type=int)
    parser_stats.add_argument('-r','--readsthreshold',dest='counts',
                              help='Min read counts',default=30,type=int)
    parser_stats.add_argument('-d','--minksdistance',dest='minks',
                              help='Minimum KS Distance for plots',default=0.0,type=float)
    parser_stats.add_argument('-p','--pvalue',dest='pvalue',
                              help='P-value for significance',default=0.01,type=float)
    parser_stats.add_argument('-o','--sort',dest='sort',
                              help='Sort data by KS of this sample')    
    
    parser_stats.set_defaults(func=run_subtool)

    args = parser.parse_args()        
    args.func(parser,args)
    
if __name__ == "__main__":
    main()
