# Pipeline for Topological Toolkit Morse-Smale Complex / Sublevelset Persistent Homology
This is a python pipeline for interfacing with TTK and run Morse-Smale Complex/Sublevelset Persistent Homology analysis

For creating VTK file you need:
  Mayavi (tvtk module): https://docs.enthought.com/mayavi/mayavi/installation.html

For running TTK you need:
  Installing TTK: https://anaconda.org/conda-forge/topologytoolkit
  Installing vtk: https://anaconda.org/conda-forge/vtk

For visualization you need:
  Matplotlib: https://anaconda.org/conda-forge/matplotlib
  For 3D iso-surface plot, you also need scikit-image (skimage module): https://scikit-image.org/docs/stable/user_guide/install.html

If there are conflicts between TTK and other packages, you might need one conda environment only for TTK and another for creating VTK files and visualization.

Please cite the following paper if you used this pipeline:

Shao-Chun Lee, Y Z, "Interpretation of autoencoder-learned collective variables using Morse–Smale complex and sublevelset persistent homology: An application on molecular trajectories", J. Chem. Phys. 160, 144104 (2024) DOI: https://doi.org/10.1063/5.0191446

And also, remember to cite TTK paper: https://topology-tool-kit.github.io/


