"""
The AGN Metric (SFErrorMetric) run/tested in this repo.
"""
import numpy as np
from astropy.stats import mad_std
from lsst.sims.maf.metrics import BaseMetric

class SFErrorMetric(BaseMetric):
    """Structure Function (SF) Uncertainty Metric. Developed on top of LogTGaps"""
    
    def __init__(self, mag, band, timesCol='observationStartMJD', allGaps=True, 
                 units='mag', maf = False, bins=np.logspace(0, np.log10(3650), 11), 
                 weight=np.full(10, 0.1), **kwargs):
        """Init function for metric class.
        
        Args:
            mag(float): The magnitude at which to evaluate this metric.
            band(str): Using observations from which filter, 
                one of 'u', 'g', 'r', 'i', 'z', 'y'.
            timesCol(str): Time column name. Defaults to "observationStartMJD".
            allGaps(bool): Whether to use all gaps (between any two pairs of 
                observations). If False, only use consecutive paris. Defaults to True.
            units(str): Unit of this metric. Defaults to "mag".
            maf(bool): Whether to use the original method from MAF to compute delta_t; 
                if False, use the vectorized method. Defaults to False.
            bins(object): An array of bin edges. Defaults to 
                "np.logspace(0, np.log10(3650), 11)" for a total of 10 bins.
            weight(object): The weight assigned to each delta_t bin for 
                deriving the final metric. Defaluts to "weight=np.full(10, 0.1)".
            
        """
        
        # Assign metric parameters to instance object
        self.timesCol = timesCol
        self.allGaps = allGaps
        self.maf = maf
        self.bins = bins
        self.weight=weight
        self.metricName = f'SFError_{mag}_{band}'
        super(SFErrorMetric, self).__init__(col=[self.timesCol, 'visitExposureTime'], 
                                            metricName=self.metricName, units=units, **kwargs)

        
    def run(self, dataSlice, slicePoint=None):
        """Code executed at each healpix pixel to compute the metric"""
        
        # remove extremely short exposures
        dataSlice = dataSlice[dataSlice['visitExposureTime'] > 5.1]
        
        # If the total number of visits < 2, mask as bad pixel
        if dataSlice.size < 2:
            return self.badval
        
        # sort data by time column
        times = np.sort(dataSlice[self.timesCol]) 
        
        # check if use all gaps (between any pairs of observations)
        if self.allGaps:
            
            # check if using original MAF method to compute time gaps
            if not self.maf:
                # use the vectorized method
                dt_matrix = times.reshape((1, times.size)) - times.reshape((times.size, 1))
                dts = dt_matrix[dt_matrix > 0].flatten().astype(np.float16)
            else:
                # use the original method --> loop through each observations
                allDiffs = []
                for i in np.arange(1,times.size,1):
                    allDiffs.append((times-np.roll(times,i))[i:])
                dts = np.concatenate(allDiffs).astype(np.float16)
        else:
            dts = np.diff(times)
        
        # bin delta_t using provided bins; if zero pair found at any delta_t bin, 
        # replace 0 with 0.01 to avoid the exploding 1/sqrt(n) term in this metric
        result, bins = np.histogram(dts, self.bins)
        new_result = np.where(result > 0, result, 0.01)
        
        # compute photometric_error^2 population variance and population mean
        # note that variance is replaced by median_absolute_deviate^2
        # mean is replaced by median in this implementation to make it robust to 
        # outliers in simulations (e.g., dcr simulations)
        err_var = dataSlice['magErr']**2
        err_var_mu = np.median(err_var)
        err_var_std = mad_std(err_var)
        
        # compute SF error
        sf_var_dt = 2*(err_var_mu + err_var_std/np.sqrt(new_result)) 
        sf_var_metric = np.sum(sf_var_dt * self.weight)
        
        return sf_var_metric