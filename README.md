# Thingness
Evolutionary exploration of cell-level morphogenesis in 2D  

![](docs/thingness_1.gif)  
Membrane-bound dot sizes show proxy _atp_ molecule distribution, the central cell dot show cell cycle progression. Ground coloring shows _aminoacids_ proxy molecule distribution in it (it doesn't diffuse inside cells yet)

## How it works
- Cell **shape** is simulated as a balance of internal pressure (based on cell volume) and the membrane's resting length
- **Chemistry** defines 40 generic molecule types and generates a reaction graph between them. Molecules inside cells are spatially-distributed along the membrane
- Besides those 40, cells contain 3 **proxy molecule** types:
    - _atp_ - acts as a generic **energy** resource, consumed or produced in chemical reactions
    - _aminoacids_ - generic **building** material, consumed for gene expression
    - _cycle_ - upon crossing the threshold, signals the initiation of **mitosis**. Expressed as a gene
- Cell **division axis** is selected akin to how the actual mitotic spindle settles into it's orientation [[1]](#ref-1) - it spins around until the highest concentration of microtubule-binding proteins on the membrane suspend it - though the spinning is replaced by getting a normal to the first component of PCA decomposition of spatial distribution of binding proteins (it's presumed uniform along the membrane right now, so only the cell shape matters)
- **Evolution** mutates genes by altering activators/repressors/product. Though it's in a very early stage overall

## Not implemented yet
- Cell -> cell and cell <-> environment molecule diffusion and signaling  
- Molecule <-> function mapping (e.g. cell adhesion, signaling, diffusion control)
- Gene <-> molecule <-> reaction relations are sketchy

## Usage
The primary notebook is _cell_shape.ipynb_  
Results are rendered to [visdom](https://github.com/facebookresearch/visdom)

## References
1. <a name="ref-1"></a>[X. Morin, Y. Bellaïche (2011). Mitotic Spindle Orientation in Asymmetric and Symmetric Cell Divisions during Animal Development](https://www.sciencedirect.com/science/article/pii/S1534580711002462)

