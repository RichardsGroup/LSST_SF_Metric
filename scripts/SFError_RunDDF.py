import pandas as pd
import numpy as np
import os, sys

# import lsst.sim.maf moduels modules
import lsst.sims.maf.db as db
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.stackers as stackers
from lsst.sims.maf.stackers import BaseStacker
import lsst.sims.maf.plots as plots
import lsst.sims.maf.metricBundles as metricBundles
from AGNMetrics import SFErrorMetric
from AGNStacker import MagErrStacker

# import convenience functions
from opsimUtils import *

# import joblib
from joblib import Parallel, delayed

# define function to run MAF on one opsim which is easily parallelziable. 
def run_mg(run, bundleDict, dbDir, outDir, metricDataPath):
    """
    Function to run pre-defined MAF metrics on one OpSim. 
    
    Args:
        run (str): The OpSim cadence run name.
        bundleDict (dict): A dictionary of MAF metrics.
        dbDir (str): The path to the OpSim databases.
        outDir (str): The path to the resultdb databases.
        metricDataPath (str): The path to the actual metric data (.npz files). 
    """
    rt = ''
    try:
        for key in bundleDict:
            bundleDict[key].setRunName(run)

        # init connection given run name
        opSimDb, resultDb = connect_dbs(dbDir, outDir, dbRuns=[run])
        # make a group
        metricGroup = metricBundles.MetricBundleGroup(bundleDict, opSimDb[run], 
                                                      metricDataPath, 
                                                      resultDb[run], verbose=False)
        metricGroup.runAll()
    
        # close sql db
        opSimDb[run].close()
        resultDb[run].close()
        
    except Exception as e:
        print(f'{run} failed!')
        print(e)
        print('----------------------')
        rt = run
            
    return rt

# function to run entire fbs version
def run_fbs(version, dbDir, outDir, metricDataPath):
    
    # create if not exists
    if not os.path.exists(os.path.abspath(outDir)):
        os.makedirs(os.path.abspath(outDir))
    
    if not os.path.exists(os.path.abspath(metricDataPath)):
        os.makedirs(os.path.abspath(metricDataPath))
            
    # get all runs
    dbRuns = show_opsims(dbDir)[:]

    # placeholder for joblib returned result
    rt = []
    rt = Parallel(n_jobs=14)(delayed(run_mg)(run, bundleDict, dbDir, outDir, metricDataPath) 
                             for run in dbRuns)

    # check failed 
    failed_runs = [x for x in rt if len(x) > 0]

    with open(f'v{version}.log', 'a') as f:
        for run in failed_runs:
            f.write(run+'\n')

#     notify.send(f"Done with FBS_v{version}!")
    
    
if __name__ == "__main__":
    
    # FBS versions to run
    versions = ['1.5', '1.6', '1.7']
    
    # get input from command line
    dbDir_temp = '/home/idies/workspace/lsst_cadence/FBS_{}/'
    outputFolder = sys.argv[1]
    
    outDir = os.path.join(outputFolder, 'ResultDBs')
    metricDataPath = os.path.join(outputFolder, 'MetricData')
    
    for version in versions:
        dbDir = dbDir_temp.format(version)
        run_fbs(version, dbDir, outDir, metricDataPath)