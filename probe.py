import numpy as np
from coefficient import ProbeBase, TMAX


class Probe(ProbeBase):
    def get_mu(self, rs, thetas):
        # TODO
        return np.ones_like(rs) * TMAX

    def get_lc(self, rs, thetas, ts):
        # TODO
        return np.ones_like(rs)
