'''
Return the value of the probe function, based on histogram.h5.
'''

import h5py
import numpy as np
from coefficient import ProbeBase

# PMT数量
N = 17612
# 液闪区域半径(mm)
R0 = 17710
# probe函数积分上限(ns)
T_MAX = 1000

class Probe(ProbeBase):
    '''
    Return the value of the probe function, based on histogram.h5.
    '''

    probe = None

    @classmethod
    def load_data(self):
        '''
        Load the probe data from the histogram.h5 file if it's not already loaded.
        This method reads the Bins and T_Bins attributes from the file and 
        stores the data in the `probe` class variable.
        '''
        if self.probe is None:
            with h5py.File('./histogram.h5', 'r') as h5file_r:
                probe_data = h5file_r["Probe"]
                # r格子数目 = θ格子数目 = self.bins
                self.bins = probe_data.attrs.get('Bins')
                # t格子数目 = self.tbins
                self.tbins = probe_data.attrs.get('T_Bins')
                # r格子
                self.r_bins = probe_data.attrs.get('R_Bins')
                # t格子数目 = self.tbins
                self.theta_bins = probe_data.attrs.get('Theta_Bins')
                self.probe = probe_data[:]

    def get_mu(self, rs, thetas):
        self.load_data()

        r_grid = np.clip(np.searchsorted(rs, self.r_bins)-1, 0, self.bins-1)
        theta_grid = np.clip(np.searchsorted(thetas, self.theta_bins)-1,
                                                            0, self.bins-1)
        return (np.sum(self.probe, axis=2))[r_grid, theta_grid] * T_MAX / self.tbins

    def get_lc(self, rs, thetas, ts):
        self.load_data()

        ts_cliped = np.clip(ts, 0, (1 - 1e-15) * T_MAX)

        r_grid = np.clip(np.searchsorted(self.r_bins, rs)-1, 0, self.bins-1)
        theta_grid = np.clip(np.searchsorted(self.theta_bins, thetas)-1,
                                                            0, self.bins-1)
        t_grid = np.floor(ts_cliped * self.tbins / T_MAX).astype(int)

        lc = self.probe[r_grid, theta_grid, t_grid]
        return lc