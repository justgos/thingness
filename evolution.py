import json
import numpy as np
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon

import chemistry
import genetics
import cytology


max_population = 1024


def sort_genomes(genomes):
    sorted_genomes = sorted(genomes, key=lambda x: x['createdTime'], reverse=True)[:max_population]
    genome_weights = np.arange(len(sorted_genomes))[::-1].astype(np.float32) + 1
    genome_weights = genome_weights / np.sum(genome_weights)
    return sorted_genomes, genome_weights


def sample_batch(epoch, env_size, ground_level, genomes, genome_map, sorted_genomes, genome_weights, genes, gene_map):
    # Setup the environment
    env = np.zeros((env_size[0], env_size[1], len(chemistry.molecule_map))).astype(np.float32)
    # Earth starts full of nutrients
    env[:, :ground_level, chemistry.molecule_map['aminoacids']] = np.expand_dims(np.linspace(1.0, 0.1, ground_level), axis=0)

    # Setup the initial cells
    batch_size = 1
    cells = []
    # Sample new genomes
    parent_batch = np.random.choice(sorted_genomes, size=batch_size, p=genome_weights)
    cell_centers = [np.array([x, ground_level]) for x in np.linspace(
        0.0, env.shape[0], batch_size, endpoint=False) + float(env.shape[0]) / batch_size * 0.5]
    for i in np.arange(batch_size):
        new_genome = parent_batch[i]
        # Mutate half of the fresh batch
#       if i < batch_size / 2:
        new_genome = genetics.mutate_genome(new_genome, genes, gene_map)
        new_genome['createdTime'] = int(epoch)
    #     print(json.dumps(new_genome))
        new_genome_id = 0
        # Check whether the mutate version already exists
        if json.dumps(new_genome['genes']) not in genome_map:
            genomes.append(new_genome)
            new_genome_id = len(genomes) - 1
            genome_map[json.dumps(new_genome['genes'])] = new_genome_id
        else:
            new_genome_id = genome_map[json.dumps(new_genome['genes'])]
        new_cell = cytology.create_cell()
        new_cell['membrane'] += cell_centers[i]
        new_cell['genome'] = new_genome_id
        new_cell['contents'][:, chemistry.molecule_map['atp']] = 300.0
        new_cell['contents'][:, chemistry.molecule_map['aminoacids']] = 100.0
        cells.append(new_cell)
    return env, cells
