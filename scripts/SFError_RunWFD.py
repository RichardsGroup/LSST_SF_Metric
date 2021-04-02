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
    
    # create a bundle dict
    bundleDict = {}
    my_bins = np.logspace(-2, np.log10(3650), 16)
    my_weights = np.full(15, 1/15)
    for gmag in [24]:

        # declare metric, slicer and sql contraint
        sf_metricG = SFErrorMetric(gmag, 'g', bins=my_bins, weight=my_weights)
        slicer = slicers.HealpixSlicer(nside=64)
        constraintG = 'filter = "g"'
        constraintG += ' and note not like "DD%"'
        constraintG += ' and proposalId = 1'

        # make a bundle
        SF_mbG = metricBundles.MetricBundle(sf_metricG, slicer, constraintG, 
                                            stackerList=[MagErrStacker(gmag)])
        summaryMetrics = [metrics.MedianMetric(), metrics.MeanMetric(), 
                          metrics.RmsMetric()]
        SF_mbG.setSummaryMetrics(summaryMetrics)

        # declare u band metric
        umag = gmag + 0.15
        sf_metricU = SFErrorMetric(umag, 'u', bins=my_bins, weight=my_weights)
        constraintU = 'filter = "u"'
        constraintU += ' and note not like "DD%"'
        constraintU += ' and proposalId = 1'

        # make a bundle
        SF_mbU = metricBundles.MetricBundle(sf_metricU, slicer, constraintU, 
                                            stackerList=[MagErrStacker(umag)])
        SF_mbU.setSummaryMetrics(summaryMetrics)

        # put into dict
        bundleDict[sf_metricG.metricName] = SF_mbG
        bundleDict[sf_metricU.metricName] = SF_mbU
        
    # get all runs
    dbRuns = show_opsims(dbDir)[:]

    # placeholder for joblib returned result
    rt = []
    rt = Parallel(n_jobs=14)(delayed(run_mg)(run, bundleDict, dbDir, outDir, metricDataPath) 
                             for run in dbRuns)
    
    # check failed 
    failed_runs = [x for x in rt if len(x) > 0]

    # rerun failed ones caused sql I/O error
    for run in failed_runs:
        print(f'Rerun failed: {run}')
        print('-------------------------------------')
        try:
            run_mg(run, bundleDict, dbDir, outDir, metricDataPath)
            failed_runs.remove(run)
        except:
            continue
            
    if len(failed_runs) > 0:
        with open(f'v{version}_SF_WFD.log', 'a') as f:
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