

def output_jgsw(rnet, hfun, tags, rsph=1.):
    """
    Output filtered reach network for mesh-gen. using jigsaw.

    """
    import numpy as np
    import jigsawpy  # here, so jigsaw is not a "strict" dep.

    geom = jigsawpy.jigsaw_msh_t()
    
    geom.mshID = "ellipsoid-mesh"
    geom.radii = \
        rsph * np.ones(3, dtype=geom.REALS_t)

    nset = []; eset = []; last = 0
    for rdat in rnet:
        if (rdat.flag >= +1):

            npts = len(rdat.vert)

            temp = jigsawpy.jigsaw_msh_t()
            temp.vert2 = np.zeros(
                (npts + 0), dtype=temp.VERT2_t)
            temp.edge2 = np.zeros(
                (npts - 1), dtype=temp.EDGE2_t)

            temp.edge2["IDtag"][:] = rdat.hpos

            indx = np.arange(0, npts - 1) + last

            last = last + npts

            temp.vert2["coord"][:, 0] = \
                [x.xlon for x in rdat.vert]

            temp.vert2["coord"][:, 1] = \
                [x.ylat for x in rdat.vert]

            temp.edge2["index"][:, 0] = indx + 0
            temp.edge2["index"][:, 1] = indx + 1

            nset.append(temp.vert2)
            eset.append(temp.edge2)
        
    geom.vert2 = np.concatenate(nset, axis=+0)
    geom.edge2 = np.concatenate(eset, axis=+0)

    geom.vert2["coord"] *= np.pi / 180.
    
    __, fmap, rmap = np.unique(
        geom.vert2["coord"], 
        return_index=True, return_inverse=True, axis=0)

    geom.vert2 = geom.vert2[fmap]
    geom.edge2["index"] = rmap[geom.edge2["index"]]

    jigsawpy.savemsh(tags + "_geom.msh", geom)
    jigsawpy.savevtk(tags + "_geo2.vtk", geom)

    geom.vert3 = np.zeros(
        geom.vert2.size , dtype=geom.VERT3_t)
    geom.vert3["coord"] = jigsawpy.S2toR3(
        geom.radii, geom.vert2["coord"])
    geom.vert2 = None

    jigsawpy.savevtk(tags + "_geo3.vtk", geom)


    init = jigsawpy.jigsaw_msh_t()

    init.mshID = "ellipsoid-mesh"
    init.radii = \
        rsph * np.ones(3, dtype=init.REALS_t)

    nset = []; eset = []; last = 0
    for rdat in rnet:
        if (rdat.flag >= 1 
                and rdat.vert[+0].seed >= 1):

            npts = +1
            temp = jigsawpy.jigsaw_msh_t()
            temp.vert2 = np.zeros(
                (npts + 0), dtype=temp.VERT2_t)
           
            temp.vert2["coord"][0, 0] = \
                rdat.vert[+0].xlon
            temp.vert2["coord"][0, 1] = \
                rdat.vert[+0].ylat

            nset.append(temp.vert2)

        if (rdat.flag >= 1 
                and rdat.vert[-1].seed >= 1):

            npts = +1
            temp = jigsawpy.jigsaw_msh_t()
            temp.vert2 = np.zeros(
                (npts + 0), dtype=temp.VERT2_t)
           
            temp.vert2["coord"][0, 0] = \
                rdat.vert[-1].xlon
            temp.vert2["coord"][0, 1] = \
                rdat.vert[-1].ylat

            nset.append(temp.vert2)
            
    init.vert2 = np.concatenate(nset, axis=+0)

    init.vert2["coord"] *= np.pi / 180.
    
    jigsawpy.savemsh(tags + "_init.msh", init)


    spac = jigsawpy.jigsaw_msh_t()

    spac.mshID = "ellipsoid-grid"
    spac.radii = \
        rsph * np.ones(3, dtype=init.REALS_t)

    spac.xgrid = hfun.xpos * np.pi / 180.
    spac.ygrid = hfun.ypos * np.pi / 180.
    spac.value = hfun.vals

    jigsawpy.savemsh(tags + "_spac.msh", spac)

    return
