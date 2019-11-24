# distutils: language=c++

import matplotlib

import os
import sys
import re
import time
import copy
import json
from collections import defaultdict
# from tqdm import tqdm
from tqdm import tqdm_notebook as tqdm

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

#from libcpp.map cimport map
#from libcpp.string cimport string
#from libcpp.vector cimport vector
#cdef map[string,int] molecule_map = chemistry.molecule_map


"""
DISCLAIMER:
This file is an attempt on converting some convert some simlation code to Cython.
Abandoned because even the 'smart' converters like "numba" will require me to
denormalize all the data almost to the level of writing GPU code, 
and it's not yet a time for that
"""


cdef dict molecule_map = chemistry.molecule_map

cpdef process_cell_reactions(float[:,:] contents, float[:,:] next_contents, int[:, :] r_reactants, int[:, :] r_products, float[:] r_rates, float[:] r_energy_costs, float sim_speed):
    # TODO: account for resource amounts becoming negative
    cdef float min_amount = 1e-6
    cdef int i
    for i in range(len(r_reactants)):
        reactants = r_reactants[i]
        products = r_products[i]
        rate = r_rates[i]
        energy_cost = r_energy_costs[i]
        # Skip reactions that do not have enough reactants
        #if np.all(np.any(contents[:, reactants] < min_amount, axis=-1)):
        #    continue
        reaction_power = np.minimum(
            np.prod(contents[:, reactants], axis=-1),
            np.min(contents[:, reactants], axis=-1)
        ) * rate
        if energy_cost > 0:
            reaction_power = np.minimum(
                reaction_power, contents[:, molecule_map['atp']] / energy_cost)
        np.add.at(
            next_contents, 
            (np.s_[:], reactants), 
            -np.expand_dims(reaction_power, axis=-1) * sim_speed)
        np.add.at(
            next_contents, 
            (np.s_[:], products), 
            np.expand_dims(reaction_power, axis=-1) * sim_speed)
        if energy_cost < 0:
            np.add.at(
                next_contents,
                (np.s_[:], molecule_map['atp']), 
                np.abs(energy_cost) * reaction_power * sim_speed)
    next_contents = np.maximum(next_contents, 0.0)
    return next_contents
