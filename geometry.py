import numpy as np


def rotate_vector(vec, angle):
    return np.array([
        [np.cos(angle), -np.sin(angle)],
        [np.sin(angle), np.cos(angle)]
    ]).dot(vec.T).T


def line_point_dist(l0, l1, p):
    return np.abs(np.cross(l1 - l0, l0 - p, axis=-1)) / np.linalg.norm(l1 - l0, axis=-1)


def get_edge_normals(m):
    edge_normals = m[(np.arange(m.shape[0]) + 1) % m.shape[0]] - m
    edge_normals = rotate_vector(edge_normals, 0.5 * np.pi)
    edge_normals /= (np.linalg.norm(edge_normals,
                                    axis=1, keepdims=True) + 1e-6)
    return edge_normals


def get_vert_normals(edge_normals):
    vert_normals = - \
        (edge_normals + edge_normals[np.arange(edge_normals.shape[0]) - 1])
    vert_normals /= (np.linalg.norm(vert_normals,
                                    axis=1, keepdims=True) + 1e-6)
    return vert_normals
