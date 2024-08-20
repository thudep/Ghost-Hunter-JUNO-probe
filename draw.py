#!/usr/bin/env python3

import numpy as np
import math
from matplotlib.colors import LogNorm
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import multiprocessing.dummy as mp
from coefficient import *
from probe import Probe
from scipy import stats
import os

matplotlib.use("agg")

plt.rcParams.update(
    {
        "font.family": "serif",  # use serif/main font for text elements
        "font.serif": ["Times"],
        "text.usetex": True,  # use inline math for ticks
        "pgf.rcfonts": False,  # don't setup fonts from rc parameters
        "font.size": 13,
    }
)

hist_rths = np.array(
    [
        (0, 0, "0"),
        (0.99, math.pi / 4, "$\\frac{\\pi}{4}$"),
        (0.99, 0, "0"),
        (0.99, math.pi, "$\\pi$"),
        (0.5, math.pi, "$\\pi$"),
        (0.7, 5 * math.pi / 6, "$\\frac 56 \\pi$"),
        (0.9, math.pi / 2, "$\\frac{\\pi}{2}$"),
    ],
    dtype=[
        ("r", np.float64),
        ("theta", np.float64),
        ("theta_text", object),
    ],
)

neighborhood_r = 0.05


def draw_pie(probe: ProbeBase, fig, ax):
    """Draw the Probe pie.

    We split the bins uniformly, and take the medium value from the bins,
    for the input :math:`r` and :math:`\\theta`.

    """
    r_bins = np.linspace(0, 1, 51)
    theta_bins = np.linspace(0, 2 * np.pi, 201)
    r_mid = (r_bins[1:] + r_bins[:-1]) / 2
    theta_mid = (theta_bins[1:] + theta_bins[:-1]) / 2
    ss = probe.get_pie(r_mid, theta_mid)
    m = ax.pcolormesh(theta_mid, r_mid, ss, norm=LogNorm(), cmap="jet")
    fig.colorbar(m)


def draw_neighborhood(fig, ax):
    """Draw the neighborhoods.

    Circle them to know where do we draw the time histograms.
    """
    mr = ax.bbox.width / 2 * neighborhood_r
    ax.scatter(
        hist_rths["theta"],
        hist_rths["r"],
        facecolors="none",
        edgecolors="k",
        s=mr**2,
        marker="o",
    )


def draw_time_hist(probe: ProbeBase, c: ConcatInfo, r, theta, fig, ax):
    """Draw :math:`R(t)`, together with the real histogram.

    .. math:: R(t)=\\frac{1}{n}\sum_{v \in \\textrm{nearby vertices}}R(r_v,\\theta_v,t)

    Where :math:`n` is the count of nearby vertices.

    The real histogram is a histogram of time of all nearby vertices-induced PEs,
    divided by :math:`n`.
    """
    sts = c.pe_ts[
        c.pe_rs**2 + r**2 - 2 * c.pe_rs * r * np.cos(c.pe_thetas - theta)
        <= neighborhood_r**2
    ]
    svf = (
        c.v_rs**2 + r**2 - 2 * c.v_rs * r * np.cos(c.v_thetas - theta)
        <= neighborhood_r**2
    )
    ts = np.linspace(0, 1000, num=10001)
    n_ts = len(ts)
    ss = probe.get_lc(np.repeat(r, n_ts), np.repeat(theta, n_ts), ts)
    ax.plot(ts, ss, label="R(t)")
    n_v = np.count_nonzero(svf)
    if n_v != 0:
        time_range = (0, 1000)
        bins = 100
        ax.hist(
            sts,
            range=time_range,
            bins=bins,
            weights=np.repeat(1.0 / n_v, len(sts)),
            label="histogram",
        )


def verf(probe: ProbeBase, c: ConcatInfo, fig, ax):
    """Draw the quotient.

    The quotient is defined as real/probe.

    See Also
    --------
    draw_pie : Calculate the probe pie
    real_pie : Calculate the real pie
    """
    N_r = int(51)
    N_θ = int(201)
    binr = np.linspace(0, 1, N_r)
    binθ = np.linspace(0, 2 * np.pi, N_θ)
    Binning = [binr, binθ]

    # to return the index of each element in 2-d histogram
    ret = stats.binned_statistic_2d(
        c.f_v_rs, c.f_v_thetas, None, "count", bins=Binning, expand_binnumbers=True
    )
    r_idx, t_idx = ret.binnumber
    # initial amplititude
    Amp = probe.get_mu(c.f_v_rs, c.f_v_thetas)
    # sum up the amplititude with the 2-d indices
    Amplitude, _, _ = np.histogram2d(
        r_idx,
        t_idx,
        bins=(np.arange(N_r) + 0.5, np.arange(N_θ) + 0.5),
        weights=Amp,
    )

    hist_PE, _, _ = np.histogram2d(c.f_pe_rs, c.f_pe_thetas, bins=Binning)

    X, Y = np.meshgrid(binθ, binr)
    cm = ax.pcolormesh(X, Y, hist_PE / Amplitude, norm=LogNorm(), cmap="jet")
    fig.colorbar(cm)


def real_pie(c: ConcatInfo, fig, ax):
    """Draw the real pie

    The pe histogram is renormalized by vertex histogram.
    """
    r_bins = np.linspace(0, 1, 51)
    theta_bins = np.linspace(0, 2 * np.pi, 201)
    Binning = [r_bins, theta_bins]
    hist_PE, binr, binθ = np.histogram2d(c.f_pe_rs, c.f_pe_thetas, bins=Binning)
    hist_predict, _, _ = np.histogram2d(c.f_v_rs, c.f_v_thetas, bins=Binning)

    X, Y = np.meshgrid(binθ, binr)
    cm = ax.pcolormesh(X, Y, hist_PE / hist_predict, norm=LogNorm(), cmap="jet")
    fig.colorbar(cm)


def Validate(probe: ProbeBase, c: ConcatInfo):
    """Calculate if the Probe is valid.

    See Also
    --------
    coefficient.ProbeBase.validate
    coefficient.ConcatInfo
    """
    return probe.validate(c.v_rs, c.v_thetas, c.pe_rs, c.pe_thetas, c.pe_ts)


def get_probe() -> ProbeBase:
    """Detemine the right probe from the coefficients.

    See Also
    --------
    coefficient.ProbeBase
    """
    return Probe()


if __name__ == "__main__":
    import argparse

    psr = argparse.ArgumentParser()
    psr.add_argument("command", type=str, help="command")
    psr.add_argument("--concat", dest="concat", type=str, help="concat file")
    psr.add_argument("-o", "--output", dest="opt", type=str, help="output file")
    args = psr.parse_args()

    if args.command == "draw":
        concat = ConcatInfo(args.concat)

        def draw_log_pie_fig(probe: ProbeBase):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection="polar", theta_offset=math.pi / 2)
            draw_pie(probe, fig, ax)
            ax.set_title(f"Pie")
            print("Pie done")
            return fig

        def draw_real_pie_fig(concat: ConcatInfo):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection="polar", theta_offset=math.pi / 2)
            real_pie(concat, fig, ax)
            draw_neighborhood(fig, ax)
            ax.set_title(f"RealPie")
            print("RealPie done")
            return fig

        def draw_quotient_fig(probe: ProbeBase, concat: ConcatInfo):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1, projection="polar", theta_offset=math.pi / 2)
            verf(probe, concat, fig, ax)
            ax.set_title(f"Quotient")
            print("Quotient done")
            return fig

        def draw_time_fig(
            probe: ProbeBase, concat: ConcatInfo, i, r, theta, theta_text
        ):
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            draw_time_hist(probe, concat, r, theta, fig, ax)
            ax.legend()
            ax.set_xlabel("t/ns")
            ax.set_title(f"Time r={r} $\\theta$={theta_text}")
            print("Time", i + 1, "done")
            return fig

        with PdfPages(args.opt) as pp:
            probe = get_probe()

            pool = mp.Pool(3 + len(hist_rths))
            figures = []

            figures.append(pool.apply_async(draw_log_pie_fig, (probe,)))
            figures.append(pool.apply_async(draw_real_pie_fig, (concat,)))
            figures.append(pool.apply_async(draw_quotient_fig, (probe, concat)))
            for i, (r, theta, theta_text) in enumerate(hist_rths):
                figures.append(
                    pool.apply_async(
                        draw_time_fig, (probe, concat, i, r, theta, theta_text)
                    )
                )

            pool.close()
            pool.join()

            for fig in figures:
                pp.savefig(figure=fig.get())

    elif args.command == "validate":
        concat = ConcatInfo(args.concat)
        probe = get_probe()
        s = Validate(probe, concat)
        if "JUNOPROBE_SCORE" in os.environ:
            t = time.time()
            print(f"{s},{t}")
            with open(os.environ["JUNOPROBE_SCORE"], mode="a") as score:
                score.write(f"{s},{t}\n")
        else:
            print(s)
    else:
        raise argparse.ArgumentError(args.command, "Invalid command")
