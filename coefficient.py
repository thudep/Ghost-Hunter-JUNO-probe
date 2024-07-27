import numpy as np
import h5py as h5
from abc import ABCMeta, abstractmethod


class ConcatInfo:
    """The loader of the concat file."""

    def __init__(self, filename: str):
        """Init ConcatInfo from an HDF5 file.

        The concat file concats vertex and PE information.
        It is useful for most calculations.

        Parameters
        ----------
        filename : str
            The file path of the concat file.
        r_max : float
            An additional and historical information, to renormalize the :math:`r`.
            It represents the real radius of the detector, and is used to ensure all the
            :math:`r` we use fall in [0, 1].
        """
        with h5.File(filename, "r", swmr=True) as file:
            concat = file["Concat"][()]
            self.pe_rs = concat["r"]
            self.pe_thetas = concat["theta"]
            self.pe_ts = concat["t"]

            self.f_pe_rs = np.hstack([self.pe_rs, self.pe_rs])
            self.f_pe_thetas = np.hstack([self.pe_thetas, 2 * np.pi - self.pe_thetas])
            self.f_pe_ts = np.hstack([self.pe_ts, self.pe_ts])

            vertices = file["Vertices"][()]
            self.v_rs = vertices["r"]
            self.v_thetas = vertices["theta"]
            self.f_v_rs = np.hstack([self.v_rs, self.v_rs])
            self.f_v_thetas = np.hstack([self.v_thetas, 2 * np.pi - self.v_thetas])


class ProbeBase(metaclass=ABCMeta):
    """The base class of all types of Probe.

    It is an abstract class, and defines the drawing and validation algorithms.
    In principle, the algorithms shouldn't be changed when a new type of Probe is added.
    Instead, just define a new class inherits `ProbeBase` and overrides the abstract methods.

    If there's really a feature request, this class should be carefully redesigned.
    """

    @abstractmethod
    def get_mu(self, rs, thetas):
        """Calculate :math:`R(r,\\theta)`

        Parameters
        ----------
        rs : numpy.ndarray
            falls in [0, 1].
        thetas : numpy.ndarray
            The same shape as `rs`. falls in [0, :math:`\\pi`]

        Returns
        -------
        numpy.ndarray
            The same shape as `rs`.
            The dimension of time is integrated.
        """
        pass

    @abstractmethod
    def get_lc(self, rs, thetas, ts):
        """Calculate :math:`R(r,\\theta,t)`.

        Parameters
        ----------
        rs : numpy.ndarray
            Falls in [0, 1].
        thetas : numpy.ndarray
            The same shape as `rs`. Falls in [0, :math:`\\pi`]
        ts : numpy.ndarray
            The same shape as `rs`. Falls in [-1, 1].

        Returns
        -------
        numpy.ndarray
            The same shape as `rs`.
        """
        pass

    def get_pie(self, rs, thetas):
        """Generate Pie from giving :math:`r` and :math:`\\theta`.

        The `rs` and `thetas` may be generated from `numpy.linspace`.

        Parameters
        ----------
        rs : numpy.ndarray
            falls in [0, 1].
        thetas : numpy.ndarray
            falls in [0, :math:`\\pi`]

        Returns
        -------
        numpy.ndarray
            The shape is `(len(thetas), len(rs))`.
        """
        thetas2, rs2 = np.meshgrid(thetas, rs)
        return self.get_mu(rs2, thetas2)

    def validate(self, v_rs, v_thetas, pe_rs, pe_thetas, pe_ts):
        """Score a Probe.

        Parameters
        ----------
        v_rs : numpy.ndarray
            The :math:`r` of vertices. falls in [0, 1].
        v_thetas : numpy.ndarray
            The :math:`\\theta` of vertices.
        pe_rs : numpy.ndarray
            The corresponding :math:`r` of PEs. falls in [0, 1].
        pe_thetas : numpy.ndarray
            The corresponding :math:`\\theta` of PEs.
        pe_ts : numpy.ndarray
            The time of PEs. falls in [0, 1].

        Returns
        -------
        float
            score.
        """
        NT = 10000
        assert np.allclose(  # check if `get_mu` and `get_lc` are consistent
            self.get_mu(v_rs, v_thetas),
            np.sum(
                self.get_lc(
                    np.tile(v_rs, (NT, 1)).T,
                    np.tile(v_thetas, (NT, 1)).T,
                    np.tile(np.linspace(0, 1, NT), (len(v_rs), 1)),
                ),
                axis=1,
            )
            / NT,
        )

        nonhit = np.sum(self.get_mu(v_rs, v_thetas))

        hit = self.get_lc(pe_rs, pe_thetas, pe_ts)

        hit = np.sum(np.log(hit))
        return hit - nonhit
