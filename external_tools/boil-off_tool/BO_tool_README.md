BOIL-OFF TOOL GUIDELINES

The Boil-off tool is written using the Cython 3.0.0a11 and Python 3.9.4 versions. In particular, the tool consists of several Python and Cython files that all contribute to building and analysing the propellant tank thermal model. The Boil-off tool package must contain the following files to properly work:

• ’Import_from_ESATAN.py’
• ’BO_tool_main_body.pyx’
• ’setup.py’
• ’Boiloff_tool.py’

In addition to these files, the 'Plots_from_main_code.py' file is added in this folder as an example on how to generate plots. Also a Radiative case from ESATAN-TMS is given here as example ('Radiative_model_GEO_single.npz'), but it can be changed according to the user's needs.

Furthermore, the following Python packages are required:

pandas
tqdm
numpy
pathlib
cython
matplotlib
CoolProp

To run the Boil-off tool, the code needs to be first compiled in C (this process is also called 'Cythonizing' the code). This is done by running the following line in the terminal:

```
python setup.py build_ext –-inplace
```

----------------------------------------------------------------------------------------------------
More guidelines can be found in:

https://repository.tudelft.nl/islandora/object/uuid:46c5872f-33df-47d2-8055-2215a9dbe9a3?collection=education
