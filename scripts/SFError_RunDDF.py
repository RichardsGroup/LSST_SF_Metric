import pandas as pd
import numpy as np
import os, sys

from notify_run import Notify
notify = Notify()

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

# add the path the scripts
sys.path.insert(0, '../scripts/')

ddfFields = ['COSMOS', 'XMM-LSS', 'ELAISS1', 'ECDFS', 'EDFS']

# define function to run MAF on one opsim which is easily parallelziable. 
def run_sf_ddf(run, src_mags, dbDir, outDir, metricDataPath, **kwargs):
    """
    Function to run SFErrorMetric on one OpSim. 
    
    Args:
        run (str): The OpSim cadence run name.
        src_mags (list): At what source magnitudes to evaluate, keys are bands. 
            Defaults to {'u': [24.15], 'g': [24]}
        dbDir (str): The path to the OpSim databases.
        outDir (str): The path to the resultdb databases.
        metricDataPath (str): The path to the actual metric data (.npz files). 
    """
    
    rt = ''
    try:
        # init connection given run name
        opSimDb, resultDb = connect_dbs(dbDir, outDir, dbRuns=[run])

        # if no DDF return
        prop_info = opSimDb[run].fetchPropInfo()
        if (not 'DD' in prop_info[1]) or (len(prop_info[1]['DD']) == 0):
            print(f'No DDF for {run}')
            return rt

        # init bundleDict
        bundleDict = {}

        # shared configs
        slicer = slicers.HealpixSlicer(nside=256)
        base_constraint = 'filter = "{}"'
        summaryMetrics = [metrics.MedianMetric(), metrics.MeanMetric(), metrics.RmsMetric()]

        # loop through bands and source mags to init metricBundle
        for band in src_mags:
            mags = src_mags[band]

            for mag in mags:

                # loop through each DDF
                for ddf in ddfFields:                
                    proposalIds = ddfInfo(opSimDb[run], ddf)['proposalId']

                    # ddf constraint based on number of fields in opsim
                    if len(proposalIds) > 1:
                        ddf_constraint = base_constraint.format(band) + \
                                        f" and (proposalId = {proposalIds[0]}" + \
                                        f" or proposalId = {proposalIds[1]})"
                    else:
                        ddf_constraint = base_constraint.format(band) + \
                                        f" and proposalId = {proposalIds[0]}"

                    sf_metric = SFErrorMetric(mag, band, **kwargs)
                    sf_metric.name = sf_metric.metricName + f'_{ddf}'
                    mb = metricBundles.MetricBundle(sf_metric, slicer, ddf_constraint, 
                                                    stackerList=[MagErrStacker(mag)])

                    bundleDict[sf_metric.name] = mb

        # set runname and summary stat
        for key in bundleDict:
            bundleDict[key].setRunName(run)
            bundleDict[key].setSummaryMetrics(summaryMetrics)
    

#     try:
        # make a group
        metricGroup = metricBundles.MetricBundleGroup(bundleDict, opSimDb[run], 
                                                      metricDataPath, 
                                                      resultDb[run], verbose=False)
        metricGroup.runAll()

        # close dbs
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
    
    # define metric parameters for DDF
    src_mags = {'u': [24.15], 'g': [24]}
    my_bins = np.logspace(-2, np.log10(3650), 21)
    my_weights = np.full(20, 1/20)

    # placeholder for joblib returned result
    rt = []
    rt = Parallel(n_jobs=15)(delayed(run_sf_ddf)(run, src_mags, dbDir, 
                                                 outDir, metricDataPath,
                                                 bins=my_bins, weight=my_weights) 
                             for run in dbRuns)

    # check failed 
    failed_runs = [x for x in rt if len(x) > 0]

    with open(f'v{version}_DDF.log', 'a') as f:
        for run in failed_runs:
            f.write(run+'\n')

    notify.send(f"Done with SF_DDF FBS_v{version}!")
    
    
if __name__ == "__main__":
    
    # FBS versions to run
    versions = ['1.5', '1.6', '1.7']
    
    # get input from command line
    dbDir_temp = '/home/idies/workspace/lsst_cadence/FBS_{}/'
    outputFolder = sys.argv[1]
    
    outDir = os.path.join(outputFolder, 'ResultDBs')
    metricDataPath = os.path.join(outputFolder, 'MetricData')
    
    for version in versions[1:]:
        dbDir = dbDir_temp.format(version)
        run_fbs(version, dbDir, outDir, metricDataPath)