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
import networkx as nx
import matplotlib.pyplot as plt

import chemistry
from plt2video import fig2data, viz_video


def hex2rgb(hex):
    h = hex.lstrip('#')
    return np.array(tuple(int(h[i:i+2], 16) for i in (0, 2, 4))).astype(np.float32) / 255


def lerp(a, b, t):
    return a + (b-a) * t


def show_frame(cells, env,):
    fig, ax = plt.subplots(figsize=np.array(env.shape[:2]) / 3)

    # env_img = np.zeros((env.shape[0], env.shape[1], 3))
    # env_img[env[:, :, chemistry.molecule_map['aminoacids']] == 0] = 1.0
    env_img = lerp(hex2rgb('#587418'), hex2rgb('#a1d32d'), 
        np.minimum(env[:, :, chemistry.molecule_map['aminoacids']:chemistry.molecule_map['aminoacids']+1], 1.0))
    env_img[env[:, :, chemistry.molecule_map['aminoacids']] == 0] = 1.0
    ax.imshow(env_img.transpose((1, 0, 2)), interpolation='bilinear', origin='lower')  # , extent=(50-0.5, (env.shape[0] - 50)-0.5, 50-0.5, (env.shape[1] - 50)-0.5)

    molecule_id = chemistry.molecule_map['atp']
    cycle_molecule_id = chemistry.molecule_map['cycle']
    for c in cells:
        # Add the cell center as an average of all the contents
        # points = np.concatenate([c['membrane'][::2], [np.mean(c['membrane'], axis=0)]])
        points = c['membrane'][::2]
        cpoints = np.array([np.mean(c['membrane'], axis=0)])
        # contents = np.power(
        #     np.concatenate([
        #         np.mean(c['contents'].reshape((c['membrane'].shape[0]//2, 2, -1)), axis=1), 
        #         [np.mean(c['contents'], axis=0)]
        #     ]), 0.5)
        contents = np.power(np.mean(c['contents'].reshape((c['membrane'].shape[0]//2, 2, -1)), axis=1), 0.5)
        ccontents = np.power([np.mean(c['contents'], axis=0)], 0.5)
    #   point_sizes = np.concatenate([300 * (c['membrane_rest_lens'][np.arange(c['membrane_rest_lens'].shape[0])-1] + c['membrane_rest_lens'] / 2), [100]])
        # point_sizes = 10 * np.concatenate([
        #     1.0 * np.ones((c['membrane_rest_lens'].shape[0]//2)), 
        #     [2.0]
        # ])
        point_sizes = 10 * 1.0 * np.ones((c['membrane_rest_lens'].shape[0]//2))
        cpoint_sizes = 10 * np.array([2.0])
    #   dots = ax.scatter(points[:,0], points[:,1], c=np.array([1.0, 0.2, 0.2, 0.0]) + np.array([0.0, 0.0, 0.0, 1.0]) * np.minimum(contents[:, :1], 0.999), s=point_sizes)
        cmap = plt.cm.get_cmap('Set2')
        dots = ax.scatter(
            points[:, 0], points[:, 1], 
            c=[cmap(molecule_id)], 
            s=point_sizes * 2.0 * np.power(contents[:, molecule_id], 0.5), 
            edgecolors='none')
        cdots = ax.scatter(
            cpoints[:, 0], cpoints[:, 1], 
            c=[cmap(cycle_molecule_id)], 
            s=cpoint_sizes * 2.0 * np.power(ccontents[:, cycle_molecule_id], 0.5), 
            edgecolors='none')
    #   dots2 = ax.scatter(points[:,0], points[:,1], c=np.array([0.2, 1.0, 0.2, 0.3]), s=point_sizes * 2.0 * np.power(contents[:, 1], 0.5), edgecolors='none')
        cell_patch = PolygonPatch(Polygon(c['membrane']), alpha=.2)  # , fc=[0.9, 0.9, 0.9]
        ax.add_patch(cell_patch)
        dots.set_clip_path(cell_patch)
        cdots.set_clip_path(cell_patch)
#     ax.add_patch(PolygonPatch(Polygon(c['membrane']), alpha=.3))
    ax.set_xlim(0, env.shape[0])
    ax.set_ylim(0, env.shape[1])
    # ax.set_xlim(-4, 4.5)
    # ax.set_ylim(-4, 4.5)
    f1 = fig2data(fig)
    plt.close()
    return f1


def render_history(cell_history, env_history):
    frames = []
    print('Rendering figures..')
    for frame_num in tqdm(np.arange(len(cell_history))):
        frames.append(show_frame(cell_history[frame_num], env_history[frame_num]))
    #     plt.close()

    frames = np.array(frames)

    # Convert frames to video and send to Visdom
    viz_opts = {
        'sim_slice': dict(
            title='Sim slice',
            width=frames.shape[2]+40,
            height=frames.shape[1]+25
        )
    }
    viz_video(frames, opts=viz_opts['sim_slice'], fps=30)


# fig, ax = plt.subplots(figsize=(8, 8))
# for c in cell_history[-1]:
#   # Add the cell center as an average of all the contents
#   points = np.concatenate([c['membrane'][::2], [np.mean(c['membrane'], axis=0)]])
#   contents = np.power(np.concatenate([np.mean(c['contents'].reshape((c['membrane'].shape[0]//2, 2, -1)), axis=1), [np.mean(c['contents'], axis=0)]]), 0.5)
# #   print(contents)
# #   point_sizes = np.concatenate([300 * (c['membrane_rest_lens'][np.arange(c['membrane_rest_lens'].shape[0])-1] + c['membrane_rest_lens'] / 2), [100]])
#   point_sizes = 300 * np.concatenate([1.0 * np.ones((c['membrane_rest_lens'].shape[0]//2)), [1.0]])
# #   dots = ax.scatter(points[:,0], points[:,1], c=np.array([1.0, 0.2, 0.2, 0.0]) + np.array([0.0, 0.0, 0.0, 1.0]) * np.minimum(contents[:, :1], 0.999), s=point_sizes)
#   cmap = plt.cm.get_cmap('Set2')
#   dots = ax.scatter(points[:,0], points[:,1], c=cmap(0), s=point_sizes * 2.0 * np.power(contents[:, 0], 0.5), edgecolors='none')
# #   dots2 = ax.scatter(points[:,0], points[:,1], c=np.array([0.2, 1.0, 0.2, 0.3]), s=point_sizes * 2.0 * np.power(contents[:, 1], 0.5), edgecolors='none')
#   cell_patch = PolygonPatch(Polygon(c['membrane']), alpha=.2)  # , fc=[0.9, 0.9, 0.9]
#   ax.add_patch(cell_patch)
#   dots.set_clip_path(cell_patch)
# #   dots2.set_clip_path(cell_patch)

#   nearby = ax.scatter(all_membranes[nearby_membranes[:,0], 0], all_membranes[nearby_membranes[:,0], 1], c=[1.0, 0.0, 0.0, 1.0], s=100, edgecolors='none')
#   nearby.set_clip_path(cell_patch)
#   nearby2 = ax.scatter(all_membranes[nearby_membranes[:,1], 0], all_membranes[nearby_membranes[:,1], 1], c=[1.0, 0.0, 0.0, 1.0], s=100, edgecolors='none')
#   nearby2.set_clip_path(cell_patch)
# ax.set_xlim(-4, 4.5)
# ax.set_ylim(-4, 4.5)
