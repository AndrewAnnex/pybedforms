## pybedforms

'pybedforms' is a Python version of 'Bedforms 4.0', a Matlab program written by David Rubin and Carissa Carter (see [this USGS publication](https://pubs.usgs.gov/of/2005/1272/) and [this interactive webpage](https://t.co/hHqVYBnGpF?amp=1)). The bedform topographies are created in exactly the same way as they were in Bedforms 4.0; however, the 3D visualization is different. We use [Mayavi](https://docs.enthought.com/mayavi/mayavi/) to build the block diagrams; these plotting functions were simplified from the ['blockdiagram' Python package](https://github.com/zsylvester/blockdiagram). All the key classes and functions are in the ['dune_topo.py' module](https://github.com/zsylvester/pybedfroms/blob/master/src/pybedforms/dune_topo.py). The ['Notebook_with_examples' jupyter notebook](https://github.com/zsylvester/pybedfroms/blob/master/Notebook_with_examples.ipynb) illustrates how to build bedform models using (1) the default parameters and (2) a set of parameters predefined by Rubin & Carter. The model can be visualized as a block diagram using the "plot_3D" function:

<img src="https://github.com/zsylvester/pybedfroms/blob/master/pybedforms_block_diagram.png" width="500">

You can also extract a "core" from the model and visualize it in 3D:

<img src="https://github.com/zsylvester/pybedfroms/blob/master/pybedforms_core.png" width="100">

## Requirements

- matplotlib
- numpy
- mayavi
- skimage
- scipy
- tqdm

