import numpy as np

molecule_types = [
    'atp',
    'aminoacids',
    'cycle',  # Proxy for cell division initiation
]
n_base_molecules = len(molecule_types)

n_aux_molecules = 20
molecule_energy = np.concatenate([
    np.ones(len(molecule_types)),
    np.clip(np.random.randn(n_aux_molecules) * 30 + 100, 10, 200),
])
molecule_types.extend(['aux{0}'.format(i) for i in np.arange(n_aux_molecules)])
aux_mol_ids = np.arange(n_base_molecules, len(molecule_types))
# _ = plt.hist(molecule_energy, bins=10)

molecule_map = {m: i for i, m in enumerate(molecule_types)}

resource_cost = {
  molecule_map['cycle']: 1.0,
}
for i, mol_id in enumerate(aux_mol_ids):
  resource_cost[mol_id] = 0.01 / (i + 1)


reactions = []
n_reactions = int(np.power(n_aux_molecules, 1.2))
for i in np.arange(n_reactions):
    n_reactants = np.clip(int(np.abs(np.random.randn()) * 0.5) + 1, 1, 2)
    reactants = np.random.choice(aux_mol_ids, n_reactants, replace=False)
    reactant_energy = np.sum([molecule_energy[m] for m in reactants])
    n_products = np.clip(int(np.abs(np.random.randn()) * 0.5) + 1, 1, 2)
    products = np.random.choice(
        [m for m in aux_mol_ids if m not in reactants.tolist()], n_reactants, replace=False)
    product_energy = np.sum([molecule_energy[m] for m in products])
    # Energy difference is either ATP consumed or produced
    unaccounted_energy = product_energy - reactant_energy
    reaction = {
        'reactants': reactants,
        'products': products,
        'rate': np.clip(np.random.power(0.05), 1e-6, 1.0),
        'energy_cost': unaccounted_energy,
    }
    reactions.append(reaction)
