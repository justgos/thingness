import matplotlib

import os
import sys
import re
import time
import copy
import json
from collections import defaultdict
from tqdm import tqdm

import numpy as np
import pandas as pd
import scipy.sparse
import scipy.spatial
from scipy.interpolate import interp1d
from shapely import affinity
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon
from descartes import PolygonPatch
from sklearn.decomposition import PCA
import networkx as nx

import geometry
import chemistry
import cytology


adhesion_max_dist = 0.2


def find_nearby_membranes(all_membranes, all_membrane_map, vert_normals):
    # Find nearby membranes for cell-cell interactions
    membrane_tree = scipy.spatial.KDTree(all_membranes)
    nearby_membranes = np.array(list(membrane_tree.query_pairs(adhesion_max_dist, p=2)))
    nearby_membrane_map = defaultdict(list)
    if nearby_membranes.shape[0] > 0:
        # Exclude same-cell membrane interactions and same-direction-facing segments
        all_vert_normals = np.concatenate(vert_normals, axis=0)
        subset = np.where(
            (all_membrane_map[nearby_membranes[:, 0], 0] !=
             all_membrane_map[nearby_membranes[:, 1], 0])
            & (np.einsum('ij,ik->i', all_vert_normals[nearby_membranes[:, 0]], all_vert_normals[nearby_membranes[:, 1]]) < 0.0)
        )
        nearby_membranes = nearby_membranes[subset]
        # {cell idx: (vert idx, other cell idx, other vert idx, 'all_membranes' vert idx)}
        for nm in nearby_membranes:
            m0 = all_membrane_map[nm[0]]
            m1 = all_membrane_map[nm[1]]
            nearby_membrane_map[m0[0]].append((m0[1], m1[0], m1[1], nm[1]))
            nearby_membrane_map[m1[0]].append((m1[1], m0[0], m0[1], nm[0]))
        nearby_membrane_map = {k: np.array(v)
                               for k, v in nearby_membrane_map.items()}
#       print(nearby_membrane_map)
    return nearby_membranes, nearby_membrane_map


def update_cell_membranes(cells):
    membrane_polys = [Polygon(cell['membrane']).buffer(0) for cell in cells]
    membrane_bounds = [p.bounds for p in membrane_polys]

    # Get normal vectors for membrane vertices
    vert_normals = [geometry.get_vert_normals(
        geometry.get_edge_normals(cell['membrane'])) for cell in cells]

    all_membranes = np.concatenate([cell['membrane'] for cell in cells], axis=0)
    # [(cell idx, vert idx), ...]
    all_membrane_map = np.concatenate([
        np.stack([
            np.repeat([i], cell['membrane'].shape[0]),
            np.arange(cell['membrane'].shape[0],)
        ], axis=1)
        for i, cell in enumerate(cells)
    ], axis=0).astype(np.int32)

    # Find inter-cell membrane vertices that are close enough for adhesion/diffusion
    nearby_membranes, nearby_membrane_map = find_nearby_membranes(
        all_membranes, all_membrane_map, vert_normals)

    # Change membrane rest length according with the cell volume
    membrane_rdists = []
    for i, cell in enumerate(cells):
        # Get all the pairwise distances between membrane vertices
        membrane_dists = scipy.spatial.distance.squareform(
            scipy.spatial.distance.pdist(cell['membrane']))
        membrane_rdists_i = 1.0 / (membrane_dists + 1e-6)
        membrane_rdists_i[np.where(membrane_dists == 0)] = 0
        membrane_rdists.append(membrane_rdists_i)

    return membrane_bounds, membrane_polys, vert_normals, \
            all_membranes, all_membrane_map, \
            nearby_membranes, nearby_membrane_map, \
            membrane_rdists


def process_cell_reactions(cell, next_cell, sim_speed):
    # TODO: account for resource amounts becoming negative
    for reaction in chemistry.reactions:
        reaction_power = np.minimum(
            np.prod(cell['contents'][:, reaction['reactants']], axis=-1),
            np.min(cell['contents'][:, reaction['reactants']], axis=-1)
        )
        if reaction['energy_cost'] > 0:
            reaction_power = np.minimum(
                reaction_power, cell['contents'][:, chemistry.molecule_map['atp']] / reaction['energy_cost'])
        np.add.at(
            next_cell['contents'], 
            (np.s_[:], reaction['reactants']), 
            -np.expand_dims(reaction_power, axis=-1) * sim_speed)
        np.add.at(
            next_cell['contents'], 
            (np.s_[:], reaction['products']), 
            np.expand_dims(reaction_power, axis=-1) * sim_speed)
        if reaction['energy_cost'] < 0:
            next_cell['contents'][:, chemistry.molecule_map['atp']
                                    ] += np.abs(reaction['energy_cost']) * reaction_power * sim_speed
    next_cell['contents'] = np.maximum(next_cell['contents'], 0.0)


def process_cell_contents_diffusion(next_cell, membrane_rdists, sim_speed):
    # TODO: non-uniform diffusion coefficients
    dcontents = np.expand_dims(
        next_cell['contents'], 1) - np.expand_dims(next_cell['contents'], 0)
#       contents_diffusion = membrane_rdists.dot(next_cell['contents']) * 0.1 * sim_speed
#       print(dcontents)
    contents_diffusion = -np.sum(np.expand_dims(membrane_rdists, -1) * dcontents, axis=1) * 0.05 * sim_speed
#       print(contents_diffusion)
    next_cell['contents'] += contents_diffusion


def process_cell_cell_collision(i, membrane, vert_normals, cell, next_cell, cells, membrane_bounds, membrane_polys, next_momentum, sim_speed):
    for j, other in enumerate(cells):
        if i == j:
            continue

        b1 = membrane_bounds[i]
        b2 = membrane_bounds[j]
        # AABB test
        if b1[2] < b2[0] or b2[2] < b1[0] or b1[3] < b2[1] or b2[3] < b1[1]:
            intersection = None
        else:
            intersection = membrane_polys[i].intersection(
                membrane_polys[j])

        # Pressure from overlapping cells
        if intersection is not None and not intersection.is_empty:
            intersects = matplotlib.path.Path(
                other['membrane']).contains_points(membrane)
            next_momentum -= vert_normals[i] * 0.5 * np.expand_dims(np.minimum(
                (intersection.area + 0.1) * intersects, 0.2 * next_cell['volume']), -1) * sim_speed


def process_membrane_tension(i, membrane, vert_normals, cell, next_cell, next_static_momentum, sim_speed):
    membrane_tension = 0.5
    to_prev = membrane[np.arange(membrane.shape[0])-1] - membrane
    to_next = membrane[(np.arange(membrane.shape[0])+1) % membrane.shape[0]] - membrane
    dp = np.linalg.norm(to_prev, axis=1)
    dn = np.linalg.norm(to_next, axis=1)
    next_static_momentum += to_prev * np.expand_dims((np.power(
        dp / next_cell['membrane_rest_lens'][np.arange(membrane.shape[0])-1], 2.0)), -1) * membrane_tension * sim_speed
    next_static_momentum += to_next * \
        np.expand_dims(
            (np.power(dn / next_cell['membrane_rest_lens'], 2.0)), -1) * membrane_tension * sim_speed
    # Internal pressure
    next_static_momentum += vert_normals[i] * np.power(
        next_cell['volume'] / (Polygon(next_cell['membrane']).area+1e-6), 2.0) * 0.1 * sim_speed


def process_cell_cell_adhesion(i, cell, all_membranes, nearby_membrane_map, adhesion_max_dist, next_momentum, sim_speed):
    if i in nearby_membrane_map:
        nearby_data = nearby_membrane_map[i]
        nearby_verts = all_membranes[nearby_data[:, 3]]
        nearby_vecs = nearby_verts - cell['membrane'][nearby_data[:, 0]]
        nearby_dists = np.linalg.norm(nearby_vecs, axis=1, keepdims=True)
        adhesion_strength = 4.0
        np.add.at(
            next_momentum, 
            nearby_data[:, 0], 
            nearby_vecs * np.power(np.maximum(adhesion_max_dist - nearby_dists, 0), 1.0) * adhesion_strength * sim_speed)


def advance_cell_division(membrane, cell, next_cell, next_static_momentum, sim_speed):
    if 'division' in cell:
        vec = membrane[cell['division']['between'][1]] - \
            membrane[cell['division']['between'][0]]
        next_static_momentum[cell['division']['between'][0]
                                ] += vec * cell['division']['stage'] * sim_speed
        next_static_momentum[cell['division']['between'][1]
                                ] -= vec * cell['division']['stage'] * sim_speed
        # FIXME: doesn't scale properly with `sim_speed`
        next_cell['division']['stage'] *= 1.0 + (0.5 * sim_speed)


def maybe_start_cell_division(next_cell):
    # Starts when the concentration of 'cycle' proxy molecule reaches the threshold
    division_start_cycle_threshold = 0.5
    if 'division' not in next_cell:
        if np.mean(next_cell['contents'][:, chemistry.molecule_map['cycle']]) >= division_start_cycle_threshold:
            cytology.start_division(next_cell)


def simulate(args):
    env, cells, steps = args
    env_history = [copy.deepcopy(env)]
    cell_history = [copy.deepcopy(cells)]
    sim_speed = 0.05
    for step in tqdm(np.arange(steps)):
        next_cells = []

        membrane_bounds, membrane_polys, vert_normals, \
            all_membranes, all_membrane_map, \
            nearby_membranes, nearby_membrane_map, \
            membrane_rdists = update_cell_membranes(cells)

        for i, cell in enumerate(cells):
            next_cell = copy.deepcopy(cell)
            next_cell['volume'] += 0.01
            membrane = cell['membrane']

            
            next_cell['membrane_rest_lens'] = np.sqrt(np.ones(
                (next_cell['membrane'].shape[0])) * next_cell['volume'] / np.pi) * 2 * np.pi / next_cell['membrane'].shape[0]

            # Process reactions inside the cell
            process_cell_reactions(cell, next_cell, sim_speed)

            # Contents diffusion
            process_cell_contents_diffusion(next_cell, membrane_rdists[i], sim_speed)

            # Cell-cell collision
            momentum_decay = 0.8
            next_momentum = np.copy(cell['membrane_momentum']) * momentum_decay
            next_static_momentum = np.zeros(next_momentum.shape)
            process_cell_cell_collision(i, membrane, vert_normals, cell, next_cell, cells, membrane_bounds, membrane_polys, next_momentum, sim_speed)

            # Membrane tension
            process_membrane_tension(i, membrane, vert_normals, cell, next_cell, next_static_momentum, sim_speed)

            # Cell-cell adhesion
            process_cell_cell_adhesion(i, cell, all_membranes, nearby_membrane_map, adhesion_max_dist, next_momentum, sim_speed)

            # Advance the division process
            advance_cell_division(membrane, cell, next_cell, next_static_momentum, sim_speed)

            # Apply the accumulated momentum changes
            next_static_momentum -= np.mean(next_static_momentum, axis=0)
            next_momentum += next_static_momentum
            next_momentum = np.clip(next_momentum, -0.1, 0.1)
            next_membrane = membrane + next_momentum

            # Division
            maybe_start_cell_division(next_cell)
            
            # Finishes when the distance between division membrane sites is less than the threshold
            next_cell['membrane'] = next_membrane
            next_cell['membrane_momentum'] = next_momentum

            division_threshold = 0.1
            if 'division' in next_cell and np.linalg.norm(next_membrane[next_cell['division']['between'][1]] - next_membrane[next_cell['division']['between'][0]]) <= division_threshold:
                c1, c2 = cytology.finish_division(next_cell)
                next_cells.extend([c1, c2])
            else:
                next_cells.append(next_cell)
        cells = next_cells
        if step % 20 == 0:
            env_history.append(copy.deepcopy(env))
            cell_history.append(copy.deepcopy(cells))

    return env_history, cell_history, all_membranes, nearby_membranes
