
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
    args.shp_file = os.path.join( # location of Hy-SHEDS
        "dat", 
        "HydroRIVERS_v10_shp", "HydroRIVERS_v10_shp",
        "HydroRIVERS_v10.shp")

    args.out_file = os.path.join( # location of filtered
        "out", 
        "HydroRIVERS_v10_shp", tag_,
        "filtered_25p0km.shp")

    args.spc_file = os.path.join( # location for spacing
        "out",
        "HydroRIVERS_v10_shp", tag_,
        "resolution.nc")

    args.msh_tags = os.path.join( # location for meshes
        "msh", 
        "HydroRIVERS_v10_msh", tag_)

    dir_ = os.path.dirname(
        os.path.abspath(args.spc_file))

    if (not os.path.exists(dir_)): os.makedirs(dir_)

    # 1. filter following jigsaw spacing
    """
    spac = jigsawpy.jigsaw_msh_t()
    name = "spac_o25_l25.msh"
    jigsawpy.loadmsh(name, spac)
    
    """

    nlon = 360; nlat = 180

    # 2. filter at const. length of 25km
    """
    data = nc.Dataset(args.spc_file, "w")
    data.createDimension("nlon", nlon)
    data.createVariable("xlon", "f8", ("nlon"))
    data["xlon"][:] = np.linspace(-180., +180., nlon)
    data.createDimension("nlat", nlat)
    data.createVariable("ylat", "f8", ("nlat"))
    data["ylat"][:] = np.linspace(-90.0, +90.0, nlat)
    data.createVariable("vals", "f4", ("nlat", "nlon"))
    data["vals"][:, :] = 25000.0  # in [m]
    data.close()
    
    """
    
    # 3. or, filter at a variable length
    data = nc.Dataset(args.spc_file, "w")
    data.createDimension("nlon", nlon)
    data.createVariable("xlon", "f8", ("nlon"))
    data["xlon"][:] = np.linspace(-180., +180., nlon)
    data.createDimension("nlat", nlat)
    data.createVariable("ylat", "f8", ("nlat"))
    data["ylat"][:] = np.linspace(-90.0, +90.0, nlat)
    data.createVariable("vals", "f4", ("nlat", "nlon"))
 
    xvec = np.asarray(data["xlon"][:]) * np.pi / 180.
    yvec = np.asarray(data["ylat"][:]) * np.pi / 180.
    
    xmat, ymat = np.meshgrid(xvec, yvec, sparse=True)
 
    # 25km globally, 5km on Mississippi basin
    filt = 25000. - 20000. * np.exp(-(
        4.0 * (xmat + 90. * np.pi / 180.) ** 2 +
        4.0 * (ymat - 30. * np.pi / 180.) ** 2) ** 4)
    
    data["vals"][:, :] = filt  # in [m]
    data.close()
    
    args.sph_size = 6371220.  # earth radius
    args.flt_endo = True  # strip no-outlet rivers??   
    args.box_xmin = -180.
    args.box_ymin = -90.0
    args.box_xmax = +180.
    args.box_ymax = +90.0

    filter_hydrosheds(args)

    # optional: build jigsaw mesh using same spacing
    """
    build_jigsaw_mesh(args)
    
    """


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



