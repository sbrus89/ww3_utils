
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import cdist


def segment(mesh):
    """
    SEGMENT: mark the cells in MESH via "connected regions",
    modifies MESH in-place such that MESH.TRIA3["IDtag"] is
    an (integer) region ID number.

    Region "boundaries" are given via the set of constraint 
    edges in the mesh: MESH.EDGE2.

    """
    # Authors: Darren Engwirda

    tidx = np.reshape(np.arange(
        0, mesh.tria3.size), (mesh.tria3.size, 1))

    ee12 = np.sort(
        mesh.tria3["index"][:, (0, 1)], axis=1)
    ee23 = np.sort(
        mesh.tria3["index"][:, (1, 2)], axis=1)
    ee31 = np.sort(
        mesh.tria3["index"][:, (2, 0)], axis=1)

    edge = np.concatenate((
        np.hstack((ee12, tidx)),
        np.hstack((ee23, tidx)),
        np.hstack((ee31, tidx))
    ))

    edge = edge[edge[:, 1].argsort(kind="stable"), :]
    edge = edge[edge[:, 0].argsort(kind="stable"), :]

    maps = np.zeros(
        (edge.shape[0] // 2, 4), dtype=edge.dtype)

    maps[:, (0, 1)] = edge[:-1:2, (0, 1)]
    maps[:, 2] = edge[:-1:2, 2]
    maps[:, 3] = edge[+1::2, 2]

    ebnd = np.sort(
        mesh.edge2["index"][:, (0, 1)], axis=1)

    bnds = np.logical_and.reduce((
        np.isin(maps[:, 0], ebnd[:, 0]),
        np.isin(maps[:, 1], ebnd[:, 1])
    ))

    maps = maps[np.logical_not(bnds.flatten()), :]

    rows = np.concatenate((maps[:, 2], maps[:, 3]))
    cols = np.concatenate((maps[:, 3], maps[:, 2]))
    vals = np.ones((rows.size), dtype=int)

    smat = csr_matrix((vals, (rows, cols)))

    cnum, cidx = connected_components(
        smat, directed=False, return_labels=True)

    mesh.tria3["IDtag"] = cidx

    # Calculate Cartesian coordinates of point in the Pacific Ocean
    R = 6371.0
    lon_pt = np.radians(-130.0)
    lat_pt = np.radians(20.0)
    x_pt = R*np.cos(lat_pt)*np.cos(lon_pt) 
    y_pt = R*np.cos(lat_pt)*np.sin(lon_pt)
    z_pt = R*np.sin(lat_pt)  
    pt = np.array([(x_pt,y_pt,z_pt)])

    # Find IDtag of element with vertex closest to ocean point
    idx = np.argmin(cdist(mesh.vert3['coord'],pt))
    elems = np.where(mesh.tria3["index"] == idx)
    ocn_id = mesh.tria3["IDtag"][elems[0][0]]

    # Create new tria3 object with elements that have the ocean IDtag
    indices, = np.where(mesh.tria3["IDtag"] == ocn_id)
    new_tria3 = []
    for el in mesh.tria3["index"][indices,:]:
      new_tria3.append((el,0))
    mesh.tria3 = np.array(new_tria3, dtype=mesh.TRIA3_t)

    return 
