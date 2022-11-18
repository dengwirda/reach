
import os
import numpy as np
import netCDF4 as nc

import jigsawpy
import argparse

from reach.filter_hydrosheds import filter_hydrosheds

class base: pass


def ex1():

#-- filter the HydroSHEDS data-set based on a nonuniform
#-- resolution function, to be used in unstructured mesh
#-- generation

    tag_ = "filtered_o25_l25"

    args = base()
    args.shp_file = os.path.join(
        "dat", 
        "HydroRIVERS_v10_shp", "HydroRIVERS_v10_shp",
        "HydroRIVERS_v10.shp")

    args.out_file = os.path.join(
        "out", 
        "HydroRIVERS_v10_shp", tag_,
        "filtered_25p0km.shp")

    args.spc_file = os.path.join(
        "out",
        "HydroRIVERS_v10_shp", tag_,
        "resolution.nc")

    args.msh_tags = os.path.join(
        "msh", 
        "HydroRIVERS_v10_msh", tag_)

    dir_ = os.path.dirname(
        os.path.abspath(args.spc_file))

    if (not os.path.exists(dir_)): os.makedirs(dir_)

    spac = jigsawpy.jigsaw_msh_t()
   #name = "spac_icom_o25_3p75_1p25_l12p5_3p75.msh"
   #name = "spac_o12p5_2p5_l6p25.msh"
    name = "spac_o25_l25.msh"
    jigsawpy.loadmsh(name, spac)

    data = nc.Dataset(args.spc_file, "w")
    data.createDimension("nlon", spac.xgrid.size)
    data.createVariable("xlon", "f8", ("nlon"))
    data["xlon"][:] = spac.xgrid * 180. / np.pi
    data.createDimension("nlat", spac.ygrid.size)
    data.createVariable("ylat", "f8", ("nlat"))
    data["ylat"][:] = spac.ygrid * 180. / np.pi
    data.createVariable("vals", "f4", ("nlat", "nlon"))
    data["vals"][:, :] = spac.value * 1000.
    data.close()

    args.sph_size = 6371220.
    args.flt_endo = True  # keep no-outlet rivers?   
    args.box_xmin = -180.
    args.box_ymin = -90.0
    args.box_xmax = +180.
    args.box_ymax = +90.0

    filter_hydrosheds(args)

    build_jigsaw_mesh(args)


def build_jigsaw_mesh(args):

#-- Based on filtered reach output, build an aligned mesh
#-- using jigsaw

    opts = jigsawpy.jigsaw_jig_t()  
    mesh = jigsawpy.jigsaw_msh_t()

    opts.geom_file = args.msh_tags + "_geom.msh"
    opts.hfun_file = args.msh_tags + "_spac.msh"
    opts.init_file = args.msh_tags + "_init.msh"
    opts.mesh_file = args.msh_tags + "_mesh.msh"
 
    opts.jcfg_file = args.msh_tags + "_opts.jig"

    opts.verbosity = +1

    opts.hfun_scal = "absolute"
    opts.hfun_hmax = float("inf")       # null HFUN limits
    opts.hfun_hmin = float(+0.00)

    opts.mesh_dims = +2
    opts.mesh_rad2 = +1.25              # relax edge-ratio
    opts.mesh_eps1 = +1.00
   #opts.mesh_top1 = True
        
    opts.optm_iter = +64                # tight optim. tol
    opts.optm_qtol = +5.E-05
    opts.optm_cost = "skew-cos"
   #opts.optm_dual = True

    jigsawpy.cmd.jigsaw(opts, mesh)

    jigsawpy.savevtk(args.msh_tags + "_mesh.vtk", mesh)


if (__name__ == "__main__"): ex1()



