import numpy as np
import matplotlib.pyplot as pl

hjd, phases, bl, sigmab, nl, sigman = np.loadtxt("grunhut_bls.txt", skiprows=1, unpack=True)
phases_fldcurve, mod_bl = np.loadtxt("fldcurv.out", skiprows=0, unpack=True)

pl.plot(phases, bl, color="black", marker=".", linestyle="none")
pl.errorbar(x=phases, y=bl, yerr=sigmab, color="black", linestyle="none")
pl.plot(phases_fldcurve, mod_bl, color="red", marker=".",
            linestyle="none")
pl.show()