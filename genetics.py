import os
import sys
import re
import time
import copy
import json

import numpy as np

import chemistry


def sample_indices(a, p):
    if len(a) > 0:
        n = np.random.binomial(len(a), p)
        if n > 0:
            return sorted(np.random.choice(np.arange(len(a)), n, replace=False), reverse=True)
    return []


def mutate_gene(gene):
    new_gene = copy.deepcopy(gene)
    add_activator_prob = 0.2
    mutate_activator_prob = 0.4
    remove_activator_prob = 0.05
    add_repressor_prob = 0.1
    mutate_repressor_prob = 0.2
    remove_repressor_prob = 0.05
    mutate_product_prob = 0.3
    # Remove/mutate/add activators
    for act_idx in sample_indices(new_gene['activators'], remove_activator_prob):
        del new_gene['activators'][act_idx]
    for act_idx in sample_indices(new_gene['activators'], mutate_activator_prob):
        new_gene['activators'][act_idx] = np.random.randint(chemistry.molecule_map['cycle'], len(chemistry.molecule_types))
    if np.random.rand() < add_activator_prob:
        new_gene['activators'].append(np.random.randint(chemistry.molecule_map['cycle'], len(chemistry.molecule_types)))
    new_gene['activators'] = sorted(new_gene['activators'])
    # Remove/mutate/add repressors
    for rep_idx in sample_indices(new_gene['repressors'], remove_repressor_prob):
        del new_gene['repressors'][rep_idx]
    for rep_idx in sample_indices(new_gene['repressors'], mutate_repressor_prob):
        new_gene['repressors'][rep_idx] = np.random.randint(chemistry.molecule_map['cycle'], len(chemistry.molecule_types))
    if np.random.rand() < add_repressor_prob:
        new_gene['repressors'].append(np.random.randint(chemistry.molecule_map['cycle'], len(chemistry.molecule_types)))
    new_gene['repressors'] = sorted(new_gene['repressors'])
    # Mutate the product
    if np.random.rand() < mutate_product_prob:
        new_gene['product'] = np.random.randint(chemistry.molecule_map['cycle'], len(chemistry.molecule_types))
    return new_gene


def mutate_genome(genome, genes, gene_map):
    new_genome = copy.deepcopy(genome)
    new_genome['fitness'] = []
    add_gene_prob = 0.2
    mutate_gene_prob = 0.2
    remove_gene_prob = 0.05
    # Remove some genes
    for gene_idx in sample_indices(new_genome['genes'], remove_gene_prob):
        del new_genome['genes'][gene_idx]
    # Mutate some genes
    for gene_idx in sample_indices(new_genome['genes'], mutate_gene_prob):
        gene_idx = np.random.randint(0, len(new_genome['genes']))
        new_gene = mutate_gene(genes[new_genome['genes'][gene_idx]])
        new_gene_id = 0
        # Check whether the mutate version already exists
        if json.dumps(new_gene) not in gene_map:
            genes.append(new_gene)
            new_gene_id = len(genes) - 1
            gene_map[json.dumps(new_gene)] = new_gene_id
        else:
            new_gene_id = gene_map[json.dumps(new_gene)]
        new_genome['genes'][gene_idx] = new_gene_id
        new_genome['genes'] = sorted(new_genome['genes'])
    # Add some genes
    if np.random.rand() < add_gene_prob:
        new_genome['genes'].append(np.random.randint(0, len(genes)))
        new_genome['genes'] = sorted(new_genome['genes'])
    return new_genome