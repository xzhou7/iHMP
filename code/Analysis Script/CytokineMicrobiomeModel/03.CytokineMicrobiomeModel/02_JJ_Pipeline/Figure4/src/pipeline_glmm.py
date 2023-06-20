"""
===========================
Jethro Johnson

Ruffus-based pipeline to run mixed-effects 
models in parallel on a PBS Torque cluster. 
===========================
"""

from ruffus import *
from ruffus.combinatorics import *

import sys
import os
import re
import time
import glob
import math
import scipy
import pickle
import gzip
import sqlite3
import collections
import shutil
import pandas as pd
import numpy as np
import random
import logging as L
import MySQLdb
import datetime
import itertools
from dateutil.relativedelta import relativedelta

import rpy2
from rpy2.robjects import r as R
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

import CGATPipelines.Pipeline as P
import CGAT.Experiment as E



# load options from the config file
PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])


###############################################################################
# Pre-process the data
###############################################################################
@transform('data/microbes.Rdata', suffix('.Rdata'), '.tsv')
def fetchMicrobeIDs(infile, outfile):
    R('''load('%(infile)s')
         write.table(df.tax, file='%(outfile)s', sep='\t', quote=FALSE)''' % locals())
    

@follows(mkdir('models.dir'))
@split(fetchMicrobeIDs,
       'models.dir/*sentinel')
def splitMicrobes(infile, outfiles):
    microbes = open(infile).readline().split()[1:]

    for microbe in microbes:
        outf = os.path.join('models.dir', microbe + '.sentinel')
        open(outf, 'w').close()


@follows(mkdir('IL17A.dir'))
@follows(mkdir('IL17F.dir'))
@follows(mkdir('IL22.dir'))
@product(splitMicrobes,
         formatter(r'models.dir/(.+).sentinel'),
         'data/IL*.Rdata',
         formatter(r'data/(.+).Rdata'),
         r'{1[1][0]}.dir/{1[0][0]}.RDS')
def runModels(infiles, outfile):
    taxon, cytokine = infiles
    taxon = P.snip(os.path.basename(taxon), '.sentinel')
    cytokine = P.snip(os.path.basename(cytokine), '.Rdata')

    df_model = infiles[1]
    df_tax = 'data/microbes.Rdata'

    log_file = P.snip(outfile, '.RDS') + '.log'
    model_file = P.snip(outfile, '.RDS') + '_brm_model.RDS'
    
    statement = ("Rscript src/run_glmm.R"
                 " %(df_model)s"
                 " %(df_tax)s"
                 " %(taxon)s"
                 " %(cytokine)s"
                 " %(outfile)s"
                 " %(model_file)s"
                 " &> %(log_file)s")
    cluster_options = "-l walltime=25:00:00,mem=10GB,nodes=1:ppn=2"
    # print(statement % locals())
    P.run()


@collate(runModels,
         regex('(.+).dir/(.+).RDS'),
         r'\1.tsv')
def collateModels(infiles, outfile):
    base = importr('base')

    results = []
    for infile in infiles:
        results.append(pandas2ri.ri2py(base.readRDS(infile)))

    df = pd.concat(results)
    df.to_csv(outfile, index=False, sep='\t')

    table_name = P.snip(outfile, '.tsv')
    df.to_sql(table_name, con=sqlite3.connect('results.sqldb'),
              if_exists='replace')

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))
