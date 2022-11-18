
import os
import time

import netCDF4 as nc
import numpy as np
import fiona
import argparse

from reach.filter import filter_init, filter_eval, \
                         filter_core, reach_dat, reach_xyz
from reach.output import output_jgsw
from reach.interp import interp_2dim, interp_grid

""" Interface to reach filtering for the HydroSHEDS
    data-set
"""
# Authors: Darren Engwirda

def filter_hydrosheds(args):

    print("Loading filter spacing...")

    ttic = time.time()

    hfun = interp_grid()
    with nc.Dataset(args.spc_file) as data:
        hfun.xpos = \
            np.asarray(data.variables["xlon"][:])
        hfun.ypos = \
            np.asarray(data.variables["ylat"][:])
        hfun.vals = \
            np.asarray(data.variables["vals"][:])

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")


    print("Loading stream network...")

    ttic = time.time()

    rnet = []; rpos = 0
    for feat in fiona.open(args.shp_file, "r"):        
        if (args.flt_endo and 
                feat["properties"]["ENDORHEIC"]):
            continue

        rdat = reach_dat()
        rdat.flag = + 1
        rdat.rpos = rpos
        rdat.rank = feat["properties"]["UPLAND_SKM"]
        rdat.hpos = feat["properties"]["HYRIV_ID"]
        rdat.dpos = feat["properties"]["NEXT_DOWN"]

        xmin, xmax = +180., -180.
        ymin, ymax = +90.0, -90.0
        for vert in feat["geometry"]["coordinates"]:
            xmin = min(xmin, vert[0])
            xmax = max(xmax, vert[0])
            ymin = min(ymin, vert[1])
            ymax = max(ymax, vert[1])

            rdat.vert.append(
                reach_xyz(vert[0], vert[1]))

        if (xmax < args.box_xmin): continue
        if (xmin > args.box_xmax): continue
        if (ymax < args.box_ymin): continue
        if (ymin > args.box_ymax): continue

        rpos = rpos + 1

        rnet.append (rdat)
        
    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")


    print("Form stream filtration...")

    ttic = time.time()

    filter_init(rnet, args.sph_size)
    filter_eval(rnet, hfun)
    filter_core(rnet, args.sph_size)

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")


    print("Writing stream network...")

    ttic = time.time()

    dir_ = os.path.dirname(os.path.abspath(args.out_file))

    if (not os.path.exists(dir_)):
        os.makedirs(dir_)

    src_ = fiona.open(args.shp_file)
    dst_ = fiona.open(
        args.out_file, "w", crs=src_.crs, schema=src_.schema, 
        driver="ESRI Shapefile")

    hmap = set()
    for rdat in rnet:
        if (rdat.flag >= +1):
            for rpos in rdat.rsub:
                hmap.add(rnet[rpos].hpos)

    for feat in src_:
        if (feat["properties"]["HYRIV_ID"] in hmap):
            dst_.write(feat)

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")

    if (args.msh_tags == ""): return


    print("Writing jigsaw outputs...")

    ttic = time.time()

    dir_ = os.path.dirname(os.path.abspath(args.msh_tags))

    if (not os.path.exists(dir_)):
        os.makedirs(dir_)

    output_jgsw(rnet, hfun, args.msh_tags, args.sph_size)

    ttoc = time.time()
    print("Done:", f"({ttoc - ttic:.2E}", "sec)\n")


