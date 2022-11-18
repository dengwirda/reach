
import random
from math import pi, sin, cos, floor

from reach.interp import interp_2dim

""" "filter" a topological river network for unstructured
    mesh generation via various proximity-type heuristics
"""
# Authors: Darren Engwirda

class reach_xyz:  # reach vertices
    __slots__ = [
        "xlon", "ylat", "xpos", "ypos", "zpos",
        "seed", "hval"]
    def __init__(self, xlon=0., ylat=0.):
    #-- [lon,lat] coordinates for vertex
        self.xlon = xlon
        self.ylat = ylat

    #-- [x, y, z] coordinates for vertex
        self.xpos = 0.
        self.ypos = 0.
        self.zpos = 0.

    #-- init. seed: >= +1 for mesh init.
        self.seed = -1

    #-- target resolution h(x) at vertex
        self.hval = 0.


class reach_dat:  # reach segments
    __slots__ = [
        "flag", "hpos", "rpos",
        "vert", 
        "xmin", "xmax", "ymin", "ymax", "zmin", "zmax",
        "rank", 
        "dpos", "upos", "rsub"]
    def __init__(self):
    #-- alive flag: >= +1 if still alive
        self.flag = -1

    #-- reach id, per external data-base
        self.hpos = -1

    #-- reach id, per internal data-sets
        self.rpos = -1

    #-- list of reach_xyz vert. in reach
        self.vert = []

    #-- [x, y, z] bounding box for reach
        self.xmin = 0.
        self.xmax = 0.
        self.ymin = 0.
        self.ymax = 0.
        self.zmin = 0.
        self.zmax = 0.

    #-- scalar priority "rank" for reach
        self.rank = 0.

    #-- rpos id for next dn-stream reach 
        self.dpos = -1
    
    #-- list of id's for up-stream reach 
        self.upos = []

    #-- list of internal merged reach id
        self.rsub = []


class reach_spt:  # simple spatial index
    __slots__ = [
        "nbox", "imax", 
        "xmin", "xmax", "ymin", "ymax", "zmin", "zmax",
        "xdel", "ydel", "zdel", 
        "ridx"]
    def __init__(self):
        self.nbox = 256  # ~50km boxes
        self.imax = 255

    #-- [x, y, z] bounding box for boxes
        self.xmin = 0.
        self.xmax = 0.
        self.ymin = 0.
        self.ymax = 0.
        self.zmin = 0.
        self.zmax = 0.

    #-- uniform [x, y, z] box resolution
        self.xdel = 0.
        self.ydel = 0.
        self.zdel = 0.

    #-- unrolled 3d list of items-in-box
        self.ridx = []


def filter_init(rnet, rsph):
    """
    Initialise reach graph: up/dn connections + x,y,z vertex
    coordinates, etc.

    """

    print("* build graphs")

    rmap = dict()
    for rdat, rpos in zip(rnet, range(len(rnet))):
        rmap[rdat.hpos] = rpos
    
    for rdat in rnet:
        rdat.rsub.append(rdat.rpos)
    
    for rdat in rnet:
        if (rdat.dpos > +0):
            rdat.dpos = rmap.get(rdat.dpos, -1)
        else:
            rdat.dpos = -1

    for rdat, rpos in zip(rnet, range(len(rnet))):
        if (rdat.dpos > -1):
            rnet[rdat.dpos].upos.append(rpos)
  
    print("* build points")

    for rdat in rnet:
        xmin, xmax = +rsph, -rsph
        ymin, ymax = +rsph, -rsph
        zmin, zmax = +rsph, -rsph
        for rpos in rdat.vert:
            xlon = rpos.xlon * pi / 180.0
            ylat = rpos.ylat * pi / 180.0
         
            rpos.xpos = rsph * cos(xlon) * \
                               cos(ylat)
            rpos.ypos = rsph * sin(xlon) * \
                               cos(ylat)
            rpos.zpos = rsph * sin(ylat)

            xmin = min (xmin, rpos.xpos)
            xmax = max (xmax, rpos.xpos)
            ymin = min (ymin, rpos.ypos)
            ymax = max (ymax, rpos.ypos)
            zmin = min (zmin, rpos.zpos)
            zmax = max (zmax, rpos.zpos)

        xmin, xmax = xmin - 1, xmax + 1
        ymin, ymax = ymin - 1, ymax + 1
        zmin, zmax = zmin - 1, zmax + 1

        rdat.xmin, rdat.xmax = xmin, xmax
        rdat.ymin, rdat.ymax = ymin, ymax
        rdat.zmin, rdat.zmax = zmin, zmax

    return


def filter_eval(rnet, hfun):
    """
    Initialise reach vert. resolution data, by interpolation
    from user-defined h(x) function.

    """

    print("* interpolator")

    xpos = [v.xlon for r in rnet for v in r.vert]
    ypos = [v.ylat for r in rnet for v in r.vert]
   
    hval = interp_2dim(hfun, xpos, ypos)
   
    hpos = 0
    for rdat in rnet:
        for vert in rdat.vert:
            vert.hval = hval[hpos]
            hpos += 1
    
    return


def filter_core(rnet, rsph):
    """
    Run core reach filtration pass: sort reaches by priority
    ranking, and delete/merge any that violate the proximity
    heuristics. 

    """

    print("* filt. stream")

    # argsort reach network by rank:
    ridx = sorted(
        range(len(rnet)), key=lambda x:rnet[x].rank)

    tree = sptree_init(rnet)

    live = set()
    for rpos in ridx: live.add(rpos)

    freq = len(ridx) / 10
    next = len(ridx) / 10

    rnum = 0
    for rpos, rnum in zip(ridx, range(len(ridx))):
        keep = True        
        rdat = rnet[rpos]

        live.discard(rpos)       

        if (rnum >= next):
            print("-", round(
                (rnum * 100.)/len(ridx), 1), "%")
            next += freq

        rnum = rnum + 1

        if (rdat.flag <= 0 or len(rdat.upos) > 1):
            continue
        
        if (keep):

    #-- 1.  Flag a reach for removal if it's "too short"...

            keep = check_short(rnet, rdat)

        if (keep):
    
    #-- 2.  Flag a reach for removal if it's "too close"...

            hmax = 0.0
            for rvrt in rdat.vert:
                hmax = max(hmax, rvrt.hval)

            halo = 1.5 * hmax
            xmin = (rdat.xmin - halo)
            xmax = (rdat.xmax + halo)
            ymin = (rdat.ymin - halo)
            ymax = (rdat.ymax + halo)
            zmin = (rdat.zmin - halo)
            zmax = (rdat.zmax + halo)

            radj = sptree_find(xmin, xmax, ymin,
                               ymax, zmin, zmax,
                               tree)

            near = list(radj & live)

            # randperm can help speed up rejections
            vert = rdat.vert
            keep = check_close(
                rnet, rdat, 
                random.sample(vert, len(vert)), 
                random.sample(near, len(near))
            )

        if (not keep):

    #-- 3.  Remove reach, and try to merge dnstream network

            for mpos in rdat.rsub:  # remove
                rnet[mpos].flag = -1
            
            dpos = rdat.dpos
            if (dpos >= +0):
                ddat = rnet[dpos]
                ddat.upos.remove(rpos)

                if (len(ddat.upos) == 1):
                    upos = ddat.upos[+0]
                    udat = rnet[upos]

                    udat.flag = +0  # merged

                    live.discard(upos)

                    xmin = min (ddat.xmin, udat.xmin)
                    xmax = max (ddat.xmax, udat.xmax)
                    ddat.xmin = xmin
                    ddat.xmax = xmax
                    
                    ymin = min (ddat.ymin, udat.ymin)
                    ymax = max (ddat.ymax, udat.ymax)
                    ddat.ymin = ymin
                    ddat.ymax = ymax

                    zmin = min (ddat.zmin, udat.zmin)
                    zmax = max (ddat.zmax, udat.zmax)
                    ddat.zmin = zmin
                    ddat.zmax = zmax

                    sptree_push(ddat.xmin, ddat.xmax,
                                ddat.ymin, ddat.ymax,
                                ddat.zmin, ddat.zmax,
                                tree, dpos)

                    ddat.vert = udat.vert[:-1] + \
                                ddat.vert

                    ddat.rsub = udat.rsub + \
                                ddat.rsub

                    ddat.upos = udat.upos

                    for npos in ddat.upos:
                        rnet[npos].dpos = dpos

    print("- 100. %")

    filter_seed(rnet, ridx, tree)

    return


def filter_seed(rnet, ridx, tree):
    """
    Flag vertices at confluences that are "well-distributed"
    wrt. the proximity heuristics. These vertices will be
    initial points in the unstructured mesh generation pass.

    """

    print("* init. vertex")

    for rpos in reversed(ridx):
        rdat = rnet[rpos]
        if (rdat.flag >= 1):

    #-- try to add up-stream vert. as a seed for mesh init.

            hrup = rdat.vert[+0].hval
            xrup = rdat.vert[+0].xpos
            yrup = rdat.vert[+0].ypos
            zrup = rdat.vert[+0].zpos

            hsqr = hrup * hrup

            hbar = hrup * 1.5
            xmin = xrup - hbar
            ymin = yrup - hbar
            zmin = zrup - hbar
            xmax = xrup + hbar
            ymax = yrup + hbar
            zmax = zrup + hbar

            radj = sptree_find(xmin, xmax, ymin, 
                               ymax, zmin, zmax,
                               tree)

            if (check_seeds(rnet, rdat, +0, list(radj))):
                rdat.vert[+0].seed = +1

    #-- try to add dn-stream vert. as a seed for mesh init.

            hrdn = rdat.vert[-1].hval
            xrdn = rdat.vert[-1].xpos
            yrdn = rdat.vert[-1].ypos
            zrdn = rdat.vert[-1].zpos

            hsqr = hrdn * hrdn

            hbar = hrup * 1.5
            xmin = xrdn - hbar
            ymin = yrdn - hbar
            zmin = zrdn - hbar
            xmax = xrdn + hbar
            ymax = yrdn + hbar
            zmax = zrdn + hbar

            radj = sptree_find(xmin, xmax, ymin, 
                               ymax, zmin, zmax,
                               tree)

            if (check_seeds(rnet, rdat, -1, list(radj))):
                rdat.vert[-1].seed = +1

    return


def stream_core(rnet, mesh):

    #assign all reaches to tree
    #loop through mesh edges and form list of reaches that intersect each edge
    
    #greedy reconstruction of flow network:
    #- find all outlet points in rnet
    #- find all mesh cells containing outlets
    #  * if cell contains multiple outlets, keep reach w largest catchment?
    #- push all outlet cells onto queue
    #- 

    return


def sptree_init(rnet):
    """
    Initialise a spatial tree based on a given reach network

    """

    print("* init. sptree")

    tree = reach_spt()
    tree.xmin = rnet[+0].xmin
    tree.xmax = rnet[+0].xmax
    tree.ymin = rnet[+0].ymin
    tree.ymax = rnet[+0].ymax
    tree.zmin = rnet[+0].zmin
    tree.zmax = rnet[+0].zmax

    for rdat, rpos in zip(rnet, range(len(rnet))):
        tree.xmin = min(tree.xmin, rdat.xmin)
        tree.xmax = max(tree.xmax, rdat.xmax)
        tree.ymin = min(tree.ymin, rdat.ymin)
        tree.ymax = max(tree.ymax, rdat.ymax)
        tree.zmin = min(tree.zmin, rdat.zmin)
        tree.zmax = max(tree.zmax, rdat.zmax)

    tree.xmin = (tree.xmin - 1)
    tree.xmax = (tree.xmax + 1)
    tree.ymin = (tree.ymin - 1)
    tree.ymax = (tree.ymax + 1)
    tree.zmin = (tree.zmin - 1)
    tree.zmax = (tree.zmax + 1)

    num = tree.nbox

    tree.xdel = (tree.xmax - tree.xmin) / num
    tree.ydel = (tree.ymax - tree.ymin) / num
    tree.zdel = (tree.zmax - tree.zmin) / num

    for lpos in range(num ** 3):
        tree.ridx.append(set())

    for rdat, rpos in zip(rnet, range(len(rnet))):
        sptree_push(rdat.xmin, 
                    rdat.xmax,
                    rdat.ymin, 
                    rdat.ymax,
                    rdat.zmin, 
                    rdat.zmax,
                    tree, rpos)

    return tree


def sptree_push(xmin, xmax, ymin, 
                ymax, zmin, zmax, 
                tree, rpos):
    """
    Push an item onto the tree based on its bounding box. If
    item is already existing, the tree is only updated where
    ths bounding box has expanded.

    """

    if (tree is None): return

    num = tree.nbox

    imin = floor((xmin - tree.xmin) / tree.xdel)
    imax = floor((xmax - tree.xmin) / tree.xdel)
    jmin = floor((ymin - tree.ymin) / tree.ydel)
    jmax = floor((ymax - tree.ymin) / tree.ydel)
    kmin = floor((zmin - tree.zmin) / tree.zdel)
    kmax = floor((zmax - tree.zmin) / tree.zdel)
    
    imax = min(tree.imax, imax)
    jmax = min(tree.imax, jmax)
    kmax = min(tree.imax, kmax)

    for ipos in range(imin, imax + 1):
        for jpos in range(jmin, jmax + 1):
            for kpos in range(kmin, kmax + 1):
                lpos = ipos * num * num + \
                       jpos * num + kpos
                tree.ridx[lpos].add(rpos)
                
    return


def sptree_find(xmin, xmax, ymin, 
                ymax, zmin, zmax, 
                tree):
    """
    Find all items that intersect with a given bounding box,
    and return as a set.

    """    

    if (tree is None): return set()

    num = tree.nbox

    imin = floor((xmin - tree.xmin) / tree.xdel)
    imax = floor((xmax - tree.xmin) / tree.xdel)
    jmin = floor((ymin - tree.ymin) / tree.ydel)
    jmax = floor((ymax - tree.ymin) / tree.ydel)
    kmin = floor((zmin - tree.zmin) / tree.zdel)
    kmax = floor((zmax - tree.zmin) / tree.zdel)
    
    imax = min(tree.imax, imax)
    jmax = min(tree.imax, jmax)
    kmax = min(tree.imax, kmax)

    ridx = set()
    for ipos in range(imin, imax + 1):
        for jpos in range(jmin, jmax + 1):
            for kpos in range(kmin, kmax + 1):
                lpos = ipos * num * num + \
                       jpos * num + kpos
                ridx.update(tree.ridx[lpos])

    return ridx


try:
#-- automagically import the compiled cython kernel routines
    from reach.kernel import check_short  # noqa
    from reach.kernel import check_close  # noqa
    from reach.kernel import check_seeds  # noqa

except ImportError:
    raise Exception("Cython kernel not found!")
