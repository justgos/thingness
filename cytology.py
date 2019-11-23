import copy

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.spatial
from scipy.interpolate import interp1d
from shapely import affinity
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon
from descartes import PolygonPatch
from sklearn.decomposition import PCA

import geometry
import chemistry


n_membrane_verts = 32
ring = np.array([(np.cos(a), np.sin(a)) for a in np.linspace(0, np.pi * 2, n_membrane_verts, endpoint=False)])


def create_cell():
    m = np.copy(ring)
    new_cell = {
        'membrane': m,
        'membrane_rest_lens': np.linalg.norm(m[(np.arange(m.shape[0]) + 1) % m.shape[0]] - m, axis=1),
        'membrane_momentum': np.zeros(m.shape),
        'membrane_elasticity': np.ones((m.shape[0])),
        'contents': np.zeros((m.shape[0], len(chemistry.molecule_types))),
        'volume': Polygon(m).area,
    }
    return new_cell


def spindle_orientation(cell):
    """Returns the maximum-variance direction of actin-binding components within the cell.
    Technically - first eigenvector from PCA of the contents-weighted, 
    normalized membrane vertex positions
    """
    m = cell['membrane'] - np.mean(cell['membrane'], axis=0)
    actin_binders_density = m / np.linalg.norm(m, axis=1, keepdims=True)
    # TODO: account for actual cell contents with actin-binding properties
    pca = PCA(n_components=2)
    pca.fit(actin_binders_density)
#   print(pca.singular_values_)
#   print(pca.components_)
    return np.array([pca.components_[0, 0], pca.components_[0, 1]])


def mitotic_cleavage_vertices(cell):
    spindle = spindle_orientation(cell)
    cleavage_vec = np.array([geometry.rotate_vector(spindle, 0.5 * np.pi)])
    membrane_cleavage_dists = geometry.line_point_dist(
        -cleavage_vec * 100.0, cleavage_vec * 100.0, cell['membrane'] - np.mean(cell['membrane'], axis=0, keepdims=True))
    cleavage_closest = np.argsort(membrane_cleavage_dists)
    # One vertex is the closes one
    div0 = cleavage_closest[0]
    # And the second is the next closest that's also somewhere on the other side of the membrane
    div1 = None
    for v in cleavage_closest[1:]:
        #     print(min(np.abs(div0 - v), np.abs(v - div0)))
        idx_dist = min(np.abs(div0 - v), np.abs(v - div0))
        if idx_dist > cell['membrane'].shape[0] * 0.25 and idx_dist < cell['membrane'].shape[0] * 0.75:
            div1 = v
            break
    return (min(div0, div1), max(div0, div1))


def start_division(cell):
    # Choose ~opposite points as the next division origins
    cell['division'] = {
        'between': mitotic_cleavage_vertices(cell),
        'stage': 0.1,
    }


def finish_division(cell):
    def make_subcell(idx):
        def tesselate(a):
            #             n = np.maximum(2*a.shape[0], 8)
            n = n_membrane_verts
#           f = interp1d(np.arange(a.shape[0]+1), np.concatenate([a, [a[-1]]]), axis=0)
            f = interp1d(np.arange(a.shape[0]), a, axis=0)
            return f(np.array(np.arange(n)) / np.float32(n) * (a.shape[0]-1))

        def cyclic(a):
            return np.concatenate([a, [a[0]]])
        c1 = copy.deepcopy(cell)
        c1['membrane'] = tesselate(cyclic(cell['membrane'])[idx])
        c1['membrane_momentum'] = tesselate(cyclic(cell['membrane_momentum'])[idx]) / 2.0
        c1['membrane_elasticity'] = tesselate(cyclic(c1['membrane_elasticity'])[idx])
        # Decrease content concentration, first x1/2 due to division, then x1/2 due to x2 tesselation
        c1['contents'] = tesselate(cyclic(c1['contents'])[idx]) / 2.0 * 0.5

        del c1['division']

        c1['volume'] = cell['volume'] / Polygon(cell['membrane']).buffer(0).area * Polygon(c1['membrane']).area * 0.5
        c1['membrane_rest_lens'] = np.sqrt(np.ones(
            (c1['membrane'].shape[0])) * c1['volume'] / np.pi) * 2 * np.pi / c1['membrane'].shape[0]
        return c1

    # Get subsets of membrane vertices that'll become the daughter cells
    c1idx = np.concatenate([np.arange(0, cell['division']['between'][0]+1),
                            np.arange(cell['division']['between'][1], cell['membrane'].shape[0]+1)])
    c2idx = np.arange(cell['division']['between'][0], cell['division']['between'][1]+1)
    return make_subcell(c1idx), make_subcell(c2idx)
