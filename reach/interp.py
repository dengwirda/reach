
import numpy as np

""" Interpolation routines for reach
"""
# Authors: Darren Engwirda

class interp_grid:
    def __init__(self):
        DIM1 = (0); DIM2 = (0, 0)
        self.xpos = \
            np.empty(DIM1, dtype=np.float64)
        self.ypos = \
            np.empty(DIM1, dtype=np.float64)
        self.vals = \
            np.empty(DIM2, dtype=np.float32)


def linear_2dim(xxww, xxee, yyss, yynn, 
                xpos, ypos,
                ffsw, ffse, ffnw, ffne):
    """
    Standard bilinear interpolation on aligned quadrilateral.

    """

    aane = (xxee - xpos) * (yynn - ypos)
    aanw = (xpos - xxww) * (yynn - ypos)
    aase = (xxee - xpos) * (ypos - yyss)
    aasw = (xpos - xxww) * (ypos - yyss)

    asum = (aane + aanw + aase + aasw) 

    return (ffnw * aase + ffne * aasw +
            ffsw * aane + ffse * aanw) / asum


def interp_2dim(ffun, xpos, ypos):
    """
    Standard bilinear interpolation on aligned tensor layout.

    """

    ipos = ffun.ypos.searchsorted(ypos)
    jpos = ffun.xpos.searchsorted(xpos)

    ipos = np.minimum(np.maximum(
        ipos, 0), ffun.ypos.size - 2)
    jpos = np.minimum(np.maximum(
        jpos, 0), ffun.xpos.size - 2)

    fval = linear_2dim(
        ffun.xpos[jpos], ffun.xpos[jpos + 1],
        ffun.ypos[ipos], ffun.ypos[ipos + 1],
        xpos, ypos,
        ffun.vals[ipos, jpos],
        ffun.vals[ipos, jpos + 1],
        ffun.vals[ipos + 1, jpos],
        ffun.vals[ipos + 1, jpos + 1])

    return fval.astype(np.float64)
