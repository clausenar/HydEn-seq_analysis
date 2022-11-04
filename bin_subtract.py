#!/usr/bin/env python3

import os, sys
import logging
import logging.config
import argparse

try:
    import numpy as np
    import pandas as pd
    from scipy import stats, integrate
    import matplotlib
    import matplotlib.pyplot as plt
except ModuleNotFoundError as e:
    print(f"\n Cannot import {e}, please activate the virtual environment\n")
    sys.exit()

import seaborn as sns

import glob
from shutil import copyfile
import itertools
import re
from datetime import datetime

import ipdb


logger = logging.getLogger(__name__)


if (os.uname()[1] == 'pinbot.biomedicine.gu.se'):
    bedgraph_pipeline_dir = "/data4/clausenLab_repository/programs/bedgraph-pipeline/"
    sys.path.append(bedgraph_pipeline_dir)
    from genome_composition import Compute
    import genome_composition
    from leading_lagging import LeadLag
    from get_base import FetchBase
    from base_class import Base
    from tools import help_functions, help_dictionaries
elif (os.uname()[1] == 'clausen-bioinf'):
    print('Loading local bedgraph-pipeline modules...')
    bedgraph_pipeline_dir = "/home/liam/Documents/code/bedgraph-pipeline/"
    sys.path.append(bedgraph_pipeline_dir)
    from genome_composition import Compute
    import genome_composition
    from leading_lagging import LeadLag
    from get_base import FetchBase
    from base_class import Base
    from tools import help_functions, help_dictionaries
else:
    sys.exit(["Cannot find the bedgraph pipeline libraries that I need, \nplease check the path and edit the source of this script"])

from base_class import Base

def setup_logging(default_path='logging.yaml',
              default_level=logging.DEBUG,
              env_key='LOG_CFG'):

    """Setup logging configuration"""

    from yaml import safe_load

    #print('Setting up logging...')

    path = default_path
    value = os.getenv(env_key, None)
    if value:
        path = value
        if os.path.exists(path):
            with open(path, 'rt') as f:
                config = safe_load(f.read())
        #print(path, config)
        logging.config.dictConfig(config)
    else:
        logging.basicConfig(level=default_level)


class ProcessBedgraphs():

    def __init__(self, fwfile, rvfile, organism, output_dir, binsize=200):

        self.logger = logging.getLogger('__init__')
        self._capture_package_versions()
        self.logger.debug("Is interactive plotting off?...{0}".format(os.environ.get("QT_QPA_PLATFORM")))

        #inherit methods and data from base_class
        self.fwfile = fwfile
        self.rvfile = rvfile
        self.organism = organism
        self.binsize = binsize
        self.output_dir = output_dir
        self._create_output_directories()

        self.accessory_data = Compute(self.organism)
        self.logger.debug("{0} {1} {2}".format(self.fwfile,
                                                       self.rvfile,
                                                       self.organism))

        self.chr_dict = self._get_chromosome_size_dict()
        self.logger.debug("Chromosome size dictionary - {0}".format(self.chr_dict))

        self.fwfile_df = self._read_in_bedgraph_file(filename=self.fwfile, sep_char="\t")
        self.rvfile_df = self._read_in_bedgraph_file(filename=self.rvfile, sep_char="\t")
        
        # check the organisms are right, otherwise difficult to debug if it's not
        assert len(set([i for i in self.chr_dict.keys()] + self.fwfile_df.iloc[:,0].unique().tolist())) == len(self.chr_dict.keys()), f'\n\nThere is a problem with chromosome names between the bedgraph files and the organism you have specified.\
            \nAre you sure you have the right organism ? \n organism --> {[i for i in self.chr_dict.keys()]} \n bedgraph file --> {self.fwfile_df.iloc[:,0].unique().tolist()} \n'

        self.logger.info("Counting the reads in the forward and reverse files")
        self.fwfile_count_total = self.count_reads_in_df(self.fwfile_df)
        self.rvfile_count_total = self.count_reads_in_df(self.rvfile_df)

        self.logger.debug("Binning forward file...{0}".format(self.fwfile))
        self.fwfile_binned_df = self.bin_bedgraph_data(dataframe=self.fwfile_df,chrdict=self.chr_dict, binsize=self.binsize)
        self.logger.debug("Binning reverse file...{0}".format(self.rvfile))
        self.rvfile_binned_df = self.bin_bedgraph_data(dataframe=self.rvfile_df,chrdict=self.chr_dict, binsize=self.binsize)

        self.logger.info("Normalising the forward and reverse reads to per million reads mapped")
        self.fwfile_normalised = self.normalise_reads_to_per_million_reads_mapped(self.fwfile_binned_df, self.fwfile_count_total, self.rvfile_count_total)
        self.rvfile_normalised = self.normalise_reads_to_per_million_reads_mapped(self.rvfile_binned_df, self.fwfile_count_total, self.rvfile_count_total)

        # self._return_data()


    def _capture_package_versions(self):

        self.logger.debug("Python version: {0}".format(sys.version))
        self.logger.debug("matplotlib version: {0}".format(matplotlib.__version__))
        self.logger.debug("pandas version: {0}".format(pd.__version__))
        self.logger.debug("numpy version: {0}".format(np.__version__))
        self.logger.debug("seaborn version: {0}".format(sns.__version__))

        return None


    def _count_reads_in_bedgraph(self, df):
        return df.iloc[:,-1].sum()

    def _read_in_bedgraph_file(self, filename, sep_char):
        # Reads a bedgraph file, skipping two rows and setting the correct variable types for the columns

        return pd.read_csv(filename,
                           sep=sep_char,
                           skiprows=2,
                           header=None,
                           dtype={0: str, 1: np.int64, 2: np.int64, 3: np.int64})

    def read_in_binned_file(self, filename, sep_char):
        return pd.read_csv(filename,
                           sep=sep_char,
                           header=None,
                           dtype={0: str, 1: np.int64, 2: np.int64, 3: np.int64})

    def write_df_csv(self, binned_df, outputfile):
        binned_df.to_csv(outputfile)


    def _get_chromosome_size_dict(self):

        self.logger = logging.getLogger('_chrsize_')
        self.logger.info('Creating chromosome size dictionary')
        chr_dict = {}
        for chromosome in help_dictionaries.numbered_chromosomes(self.organism).keys():
            # why was this added ?
        #    if len(chromosome[3:]) <= 2:
            chr_dict[chromosome] = self.accessory_data.size(chromosome)
        return chr_dict

    def _create_output_directories(self):
        self.logger = logging.getLogger('_osfunc_')
        self.logger.info('Creating output directories')
        directories = ['sample_divide_control_ratios',
                       'exploratory_strand_ratios',
                       'normalised_binned_bedgraphs',
                       'sample_minus_control_ratios_transposed',
                       'summary',
                       'sample_minus_control_ratios']

        for extra_dir in directories:
            if not os.path.exists(self.output_dir + extra_dir):
                print("Creating directory...{0}".format(extra_dir))
                os.makedirs(self.output_dir + extra_dir)
        return "True"

    def rename_columns(self, dataframe):
        return dataframe.rename(columns={0: "chromosome", 1 : "start", 2: "end", 3: "value"}, inplace=True)

    def rename_binned_columns(self, dataframe):
        return dataframe.rename(columns={0: "chromosome", 1 : "start", 2: "value"}, inplace=True)

    def change_column_dtype(self, dataframe):
        #dataframe[['start','end','value']] = dataframe[['start','end', 'value']].apply(pd.to_numeric)
        return dataframe[['start','end', 'value']].apply(pd.to_numeric)

    def define_bin_array(self, start, stop, step):
        return np.arange(start, stop, int(step))


    # def _return_bedgraph_filename_from_binned(self, filename):

    #     return(filename[:-15]+".bedgraph")


    def bin_bedgraph_data(self, dataframe, chrdict, binsize):

        self.logger = logging.getLogger('_mainbinfunc_')

        binned_df = pd.DataFrame()
        self.rename_columns(binned_df)
        self.rename_columns(dataframe)

        for chromosome, chrsize in chrdict.items():
            self.logger.info('Binning chromosome {0}'.format(chromosome))
            chr_subset = dataframe.loc[(dataframe.chromosome == chromosome)].copy()
            #chr_subset.fillna(0)
            chr_bins = self.define_bin_array(start=0, stop=chrsize, step=self.binsize)
            chr_subset.loc[:,'bins'] = pd.cut(x=chr_subset.loc[:,'start'],
                                        bins=chr_bins,
                                        labels=chr_bins[:-1],
                                        include_lowest=True)

            binned_df = pd.concat([binned_df,
                                   chr_subset.groupby(['chromosome','bins'],
                                                      axis=0,
                                                      as_index=False).value.sum()])
            binned_df = binned_df.fillna(0)
        
        return binned_df

    def count_reads_in_df(self, df):
        return df.iloc[:,3].sum()

    def normalise_reads_to_per_million_reads_mapped(self, df, fwfile_count_total, rvfile_count_total):

        columns = df.columns
        total_count=fwfile_count_total+rvfile_count_total
        df.iloc[:,-1] = 1000000*(df.iloc[:,-1]/total_count)

        return df

    def _remove_negative_values(self, value):

        if value < 0:
            return 0
        else:
            return value

    def _check_negative_values(self,df):

        df.iloc[:,-1] = df.iloc[:,-1].apply(lambda cell: self._remove_negative_values(cell))

        return df

    def divide_binned_values_in_last_df_column(self, sample, control):

        subtracted_series = sample.iloc[:,:-1].copy(deep=True)
        subtracted_series.loc[:,'fraction'] = pd.Series.divide(sample.iloc[:,-1].astype(float),
                                                           control.iloc[:,-1].astype(float))

        return subtracted_series

    def subtract_binned_values_in_last_df_column(self, sample, control):

        subtracted_df = sample.iloc[:,:-1]

        subtracted_df.loc[:,'fraction'] = pd.Series.subtract(sample.iloc[:,-1].astype(float),
                                                             control.iloc[:,-1].astype(float))

        # then make all negative values 0
        subtracted_df.loc[:,'fraction'] = subtracted_df.loc[:,'fraction'].apply(lambda value: self._remove_negative_values(value))
        return subtracted_df

    def determine_ratio_of_forward_reads_relative_to_reverse(self, fwd, rev):
        
        #if (fwd == 0) and (rev != 0):
        #    return 0
        #elif (fwd != 0) and (rev == 0):
        #    return 1
        if (pd.isnull(fwd) == True) and (pd.isnull(rev) == True):
            return "NaN"
        elif (fwd == "NaN") and (rev != "NaN"):
            return 0
        elif (fwd != "NaN") and (rev == "NaN"):
            return 1
        elif (fwd == "NaN") and (rev == "NaN"):
            return "NaN"
        elif (fwd != 0) and (rev != 0):
            return np.true_divide(fwd,fwd+rev)

    def determine_forward_fraction(self, fwd, rev):
        
        if (fwd == 0) and (rev != 0):
            return 0
        elif (fwd != 0) and (rev == 0):
            return 1
        elif (pd.isnull(fwd) == True) and (pd.isnull(rev) == True):
            return "NaN"
        elif (fwd == "NaN") and (rev != "NaN"):
            return 0
        elif (fwd != "NaN") and (rev == "NaN"):
            return 1
        elif (fwd == "NaN") and (rev == "NaN"):
            return "NaN"
        elif (fwd == 0) and (rev == 0):
            return "NaN"
        elif (fwd != 0) and (rev != 0):
            return fwd/(rev+fwd)

    def calculate_ratio_of_reads(self, fwd_df, rev_df, method):

        forward_fraction_df = pd.DataFrame()
        forward_fraction_df = fwd_df.iloc[:,0:2]

        forward_fraction_df.loc[:,'fwd'] = fwd_df.iloc[:,-1]
        forward_fraction_df.loc[:,'rev'] = rev_df.iloc[:,-1]

        # calculate_replication_fork_directionality
        try:
            if method == "subtract":
                forward_fraction_df['fraction'] = forward_fraction_df.apply(lambda row: self.determine_forward_fraction(row['fwd'], row['rev']), axis=1)
            elif method == "divide":
                forward_fraction_df['fraction'] = forward_fraction_df.apply(lambda row: self.determine_ratio_of_forward_reads_relative_to_reverse(row['fwd'], row['rev']), axis=1)
        except ValueError:
            ipdb.set_trace()

        # drop the 'fwd' and 'rev' columns, we don't need them anymore
        forward_fraction_df.drop(['fwd', 'rev'], axis=1, inplace=True)

        return forward_fraction_df

    def return_chromosome_names_as_numbered_dict(self, organism):

        return help_dictionaries.numbered_chromosomes(organism)

    def transpose_df_to_chromosome_by_column(self, df):

        self.logger = logging.getLogger('__transpose__')

        transposed_df = pd.DataFrame()
        transposed_df_description = pd.DataFrame()
        length_of_chr_series = {}

        chr_dict = self.return_chromosome_names_as_numbered_dict(self.organism)

        self.logger.debug("Transposing {0} chromosomes".format(len(chr_dict.keys())))

        for chromosome in df.iloc[:,0].unique():

            column_name = chr_dict.get(chromosome)
            series = df[df.iloc[:,0] == chromosome].iloc[:,-1].reset_index(drop=True)
            length_of_chr_series[chromosome] = len(series)

            series.name = "X" + str(column_name)

            if transposed_df.shape[0] == 0:

                transposed_df = series
                transposed_df_description = pd.DataFrame(data=series.astype(float).describe(),
                                                         index=series.astype(float).describe().index)
            else:
                transposed_df = pd.concat([transposed_df, series], axis=1)
                transposed_df_description = pd.concat([transposed_df_description, series.astype(float).describe()], axis=1)

        most_records = max(length_of_chr_series, key=length_of_chr_series.get)

        # If the sample only has one chromosome, it will only have one series
        if type(transposed_df) == pd.core.series.Series:
            transposed_df = pd.concat([transposed_df, df[df.iloc[:,0] == most_records].iloc[:,1]], axis=1).reset_index()
            transposed_df.drop(['index'], axis=1, inplace=True)
            transposed_df.rename(index=str, columns={'bins':'Position'}, inplace=True)

        else:
            transposed_df['Position'] = df[df.iloc[:,0] == most_records].iloc[:,1].reset_index(drop=True)

        column_order = ["Position"] + ["X"+str(i) for i in range(1,len(chr_dict.keys())+1)]

        transposed_df = transposed_df.reindex(column_order, axis=1)

        return transposed_df, transposed_df_description

    def get_normalised_binned_data(self):

        return self.fwfile_normalised, self.rvfile_normalised

    def get_binned_reads(self):

        return self.fwfile_binned_df, self.rvfile_binned_df

def create_output_filenames(output_dir, fwd, rev, ctrl_fwd, ctrl_rev, binsize):


    fwd = os.path.basename(fwd)
    rev = os.path.basename(rev)
    ctrl_fwd = os.path.basename(ctrl_fwd)
    ctrl_rev = os.path.basename(ctrl_rev)

    sample_bin_fwd = ''.join([output_dir,"normalised_binned_bedgraphs/",fwd[:-8],"tsv"])
    sample_bin_rev = ''.join([output_dir,"normalised_binned_bedgraphs/",rev[:-8],"tsv"])

    control_bin_fwd = ''.join([output_dir,"normalised_binned_bedgraphs/",ctrl_fwd[:-8]+"tsv"])
    control_bin_rev = ''.join([output_dir,"normalised_binned_bedgraphs/"+ctrl_rev[:-8]+"tsv"])

    exploratory_sample_ratio = ''.join([output_dir, "exploratory_strand_ratios/", fwd, "_AND_", rev, "_STRAND_RATIOS.tsv"])

    exploratory_ctrl_ratio = ''.join([output_dir, "exploratory_strand_ratios/", ctrl_fwd, "_AND_", ctrl_rev, "_STRAND_RATIOS.tsv"])

    subtraction_ratio = ''.join([output_dir, "sample_minus_control_ratios/", fwd[:-8], "_MINUS_",
                                 ctrl_fwd[:-8], "tsv"])

    division_ratio = ''.join([output_dir, "sample_divide_control_ratios/", fwd[:-17], "_DIVIDE_",
                              ctrl_fwd[:-17], ".tsv"])
                              
    div_minmaxlintrans = ''.join([output_dir, "sample_divide_control_ratios/", fwd[:-17], "_MINMAXLIN_",
                              ctrl_fwd[:-17], ".tsv"])

    division_sqrt_ratio = ''.join([output_dir, "sample_divide_control_ratios/", fwd[:-17], "_DIVIDE_SQRT_",
                              ctrl_fwd[:-17], ".tsv"])

    division_sqrt_norm_ratio = ''.join([output_dir, "sample_divide_control_ratios/", fwd[:-17], "_DIVIDE_SQRT_NORM_",
                                        ctrl_fwd[:-17], ".tsv"])

    transposed = ''.join([output_dir, "sample_minus_control_ratios_transposed/", fwd[:-17], "_STRAND_RATIO_", ctrl_fwd[:-17], "_", binsize, "bp_bin_transposed.tsv"])

    transposed_summary = ''.join([output_dir, "summary/", fwd[:-17], "_MINUS_", ctrl_fwd[:-17], "_summary.tsv"])

    transposed_div = ''.join([output_dir, "sample_divide_control_ratios_transposed/", fwd[:-17], "_STRAND_RATIO_", ctrl_fwd[:-17], "_", binsize, "bp_bin_transposed.tsv"])

    transposed_div_summary = ''.join([output_dir, "summary/", fwd[:-17], "_DIV_", ctrl_fwd[:-17], "_div_summary.tsv"])

    return (sample_bin_fwd,
            sample_bin_rev,
            control_bin_fwd,
            control_bin_rev,
            exploratory_sample_ratio,
            exploratory_ctrl_ratio,
            subtraction_ratio,
            division_ratio,
            div_minmaxlintrans,
            division_sqrt_ratio,
            division_sqrt_norm_ratio,
            transposed,
            transposed_summary,
            transposed_div_summary)

def min_max_linear_transformation(df):
    mm_df = df.copy(deep=True)
    # exclude mitochondria from normalisation
    mm_df = mm_df[mm_df.loc[:,'chromosome'] != sample.accessory_data.get_mito_name()]
    mm_df.iloc[:,-1] = (mm_df.iloc[:,-1] - mm_df.iloc[:,-1].min()) / (mm_df.iloc[:,-1].max() - mm_df.iloc[:,-1].min())
    return mm_df


def normalise_sqrt_ratio_to_between_zero_and_one(df):
    norm_df = df.copy(deep=True)

    # exclude mitochondria from normalisation
    norm_df = norm_df[norm_df.loc[:,'chromosome'] != sample.accessory_data.get_mito_name()]
    
    norm_df.iloc[:,-1] = 1/(1+(1/norm_df.iloc[:,-1]))
    return norm_df

def sqrt_ratio_column(df):
    sqrt_df = df.copy(deep=True)
    sqrt_df = sqrt_df[sqrt_df.loc[:,'chromosome'] != sample.accessory_data.get_mito_name()]
    sqrt_df.loc[:,'fraction'] = np.sqrt(sqrt_df.loc[:,'fraction'])
    return sqrt_df


if __name__ == '__main__':

    '''
        The method that runs when the script is invoked from the commandline
    '''
    import logging
    import logging.config
    import os
    parser = argparse.ArgumentParser(description='\nThis script reads in the two bedgraph files from a particular directory, bins them, and calculates the replication fork directionality')

    parser.add_argument('--p',
                        metavar='- directory path containing bedraph files',
                        required=True,
                        action='store',
                        help='Directory containing bedgraph files either linked or hard linked.')

    parser.add_argument('--1',
                        metavar='fwd_file',
                        required=True,
                        action='store',
                        help='Name of the FORWARD bedgraph file containing the SIGNAL.')

    parser.add_argument('--2',
                        metavar='rev_file',
                        required=True,
                        action='store',
                        help='Name of the REVERSE bedgraph file containing the SIGNAL.')

    parser.add_argument('--3',
                        metavar='fwd_file_ctrl',
                        required=True,
                        action='store',
                        help='Name of the FORWARD bedgraph file containing the CONTROL.')

    parser.add_argument('--4',
                        metavar='rev_file_ctrl',
                        required=True,
                        action='store',
                        help='Name of the REVERSE bedgraph file containing the CONTROL.')


    parser.add_argument('--o',
                        metavar='output directory',
                        required=True,
                        action='store')

    parser.add_argument('--a',
                        metavar='organism',
                        required=True,
                        action='store',
                        help='Choose an organism')

    parser.add_argument('--b',
                        metavar='bin size',
                        action='store',
                        help='The size of the bp window in which to bin the reads from the bedgraph files (default=200) - optional')


    try:
        results = parser.parse_args(None if sys.argv[1:] else ['-h'])

        if not (results.p or results.o or results.a):
            parser.error("You are missing some parameters.")
            print(f'Parameters found {results.keys()}') 
            parser.print_help()
        else:
        #    if check_if_path_exists(results.o) == False:
        #        print("Output path does not exist")
            setup_logging()
            os.environ["QT_QPA_PLATFORM"] = "offscreen"

            cmdline_dict = vars(results)
            output_dir = cmdline_dict.get('o')
            data_dir = cmdline_dict.get('p')
            organism = cmdline_dict.get('a')
            fwfile = os.path.join(data_dir, cmdline_dict.get('1'))
            rvfile = os.path.join(data_dir, cmdline_dict.get('2'))
            control_fwfile = os.path.join(data_dir, cmdline_dict.get('3'))
            control_rvfile = os.path.join(data_dir, cmdline_dict.get('4'))


            if len({fwfile,rvfile,control_fwfile,control_rvfile}) == 4:
                pass
            else:
                print('The filenames are not unique, please check')
                sys.exit(-1)


            if cmdline_dict.get('b') is None:
                binsize = 200
            else:
                binsize = cmdline_dict.get('b')

            print('Creating filenames')

            (sample_fwd_bin_filename,
             sample_rev_bin_filename,
             control_fwd_bin_filename,
             control_rev_bin_filename,
             sample_ratio_filename,
             ctrl_ratio_filename,
             subtraction_ratio_filename,
             division_ratio_filename,
             div_ratio_linstrans_filename,
             division_sqrt_filename,
             division_sqrt_norm_ratio_filename,
             transposed_filename,
             transposed_summary_filename,
             transposed_div_filename,
             transposed_div_summary_filename) = create_output_filenames(output_dir,
                                                                        fwfile,
                                                                        rvfile,
                                                                        control_fwfile,
                                                                        control_rvfile,
                                                                        binsize)

            print("Processing sample files")
            sample = ProcessBedgraphs(fwfile=fwfile,
                                      rvfile=rvfile,
                                      organism=organism,
                                      output_dir=output_dir,
                                      binsize=binsize)

            # comes as a tuple
            sample_data_fwd, sample_data_rev = sample.get_normalised_binned_data()

            sample_data_fwd.to_csv(sample_fwd_bin_filename, sep="\t", index=False)
            sample_data_rev.to_csv(sample_rev_bin_filename, sep="\t", index=False)

            print("Processing control files")
            control = ProcessBedgraphs(fwfile=control_fwfile,
                                       rvfile=control_rvfile,
                                       organism=organism,
                                       output_dir=output_dir,
                                       binsize=binsize)
            # comes as tuple
            control_data_fwd, control_data_rev = control.get_normalised_binned_data()

            control_data_fwd.to_csv(control_fwd_bin_filename, sep="\t", index=False)
            control_data_rev.to_csv(control_rev_bin_filename, sep="\t", index=False)




            #
            # Determine the forward strand ratio using read subtraction
            #


            subtraction_fwd = sample.subtract_binned_values_in_last_df_column(sample=sample_data_fwd,
                                                                              control=control_data_fwd)

            subtraction_rev = sample.subtract_binned_values_in_last_df_column(sample=sample_data_rev,
                                                                              control=control_data_rev)
            
            subtraction_ratios = control.calculate_ratio_of_reads(fwd_df=subtraction_fwd,
                                                                  rev_df=subtraction_rev,
                                                                  method="subtract")

            subtraction_ratios.to_csv(subtraction_ratio_filename, sep="\t", index=False)

            # Transpose the subtraction data for the R script plotting function

            transposed_data, transposed_data_description = sample.transpose_df_to_chromosome_by_column(subtraction_ratios)

            transposed_data.to_csv(transposed_filename, sep="\t", index=False, header=True)

            transposed_data_description.to_csv(transposed_summary_filename, sep="\t", header=True)

            #
            # Determine the strand bias using forward/reverse reads
            #

            sample_div_ratios = sample.calculate_ratio_of_reads(fwd_df=sample_data_fwd,
                                                                rev_df=sample_data_rev,
                                                                method="divide")

            sample_div_ratios.to_csv(sample_ratio_filename, sep="\t", index=False)

            control_div_ratios = control.calculate_ratio_of_reads(fwd_df=control_data_fwd,
                                                                  rev_df=control_data_rev,
                                                                  method="divide")

            control_div_ratios.to_csv(ctrl_ratio_filename, sep="\t", index=False)

            division_ratios = sample.divide_binned_values_in_last_df_column(sample=sample_div_ratios,
                                                                            control=control_div_ratios)

            transposed_div_data, transposed_div_data_description = sample.transpose_df_to_chromosome_by_column(division_ratios)

            transposed_div_data.to_csv(transposed_div_filename, sep="\t", index=False, header=True)

            transposed_div_data_description.to_csv(transposed_div_summary_filename, sep="\t", header=True)
            
            # Output the minmax transformation

            minmax_df = min_max_linear_transformation(division_ratios)
            
            #minmax_df.to_csv(div_ratio_linstrans_filename, sep="\t", index=False)
            
            minmax_df_T, minmax_df_description_T = sample.transpose_df_to_chromosome_by_column(minmax_df)

            minmax_df_T.to_csv(div_ratio_linstrans_filename, sep="\t", index=False)

            # Output the Sqrt transformation

            division_ratios_sqrt = sqrt_ratio_column(division_ratios)

            division_ratios_sqrt_T, division_ratios_sqrt_description = sample.transpose_df_to_chromosome_by_column(division_ratios_sqrt)
            
            division_ratios_sqrt_T.to_csv(division_sqrt_filename, sep="\t", index=False)

            # Calculate and output sqrt normalisation

            division_ratios_sqrt_norm = normalise_sqrt_ratio_to_between_zero_and_one(division_ratios_sqrt)
       
            division_ratios_sqrt_norm_T, division_ratios_sqrt_norm_description = sample.transpose_df_to_chromosome_by_column(division_ratios_sqrt_norm)

            division_ratios_sqrt_norm_T.to_csv(division_sqrt_norm_ratio_filename, sep="\t", index=False)

    except IOError as e:
        parser.error(str(e))
