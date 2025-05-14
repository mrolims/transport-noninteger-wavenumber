# transport-noninteger-wavenumber

Code repository accompanying the publication entitled "Transport mechanisms associated with non-integer wavenumbers in a nontwist map".

This project contains the code to generate and plot all the data from all figures.

## Requirements

To generate the all the data, you will need a Fortran compiler, either [gfortran](https://gcc.gnu.org/fortran/) or [ifx](https://www.intel.com/content/www/us/en/developer/articles/release-notes/fortran-compiler/2025.html). All figures are made with Python, using the [matplotlib](https://matplotlib.org) library. To install all required Python packages, please run 

    pip install -r requirements.txt


## Generating the data

All the data will be stored in a directory called `Data/`, so please remember to create it.

### Figures 1 and 2

The data of Figs. 1 and 2 are generated within the `Plots.ipynb` Jupyter notebook under the Fig. 1 and Fig. 2 cells.

### Figure 3

To generate the transmissivity data from Fig. 3, compile the `transmissivity_vs_m.f90` program as

    gfortran params_dp.f90 functions.f90 transmissivity_vs_m.f90 -fopenmp

or

    ifx params_dp.f90 functions.f90 transmissivity_vs_m.f90 -fopenmp

This program takes on seven parameters: `num_ic`, `m_ini`, `m_end`, `dm`, `y0`, `esc_y`, and `N`. Run the program as

    ./a.out num_ic m_ini m_end dm y0 esc_y N

with `num_ic = 10000`, `m_ini = -10`, `m_end = 10`, `dm = 0.01`, and `N = 1000000`. The `y0` and `esc_y` parameters have different values. For the blue curve, use `y0 = -1` and `esc_y = 1` and for the red curve, use `y0 = 1` and `esc_y = -1`.

### Figures 4, 5, and 6

To generate the elliptic and hyperbolic points of both the upper and lower period-11 stability islands, run

    python exe_elliptic_points.py
    python exe_hyperbolic_points.py

These scripts will compile the `elliptic_points.f90` and `hyperbolic_points.f90` programs and execute them with the appropriate parameter values. After that, run

    python exe_manifolds.py

This script will compile the `manifolds.f90` program and execute it to generate the stable and unstable manifolds of both upper and lower hyperbolic points. The `exe_*.py` scripts use `gfrotran` by default. You can use `ifx` instead by running

    python exe_elliptic_points.py ifx
    python exe_hyperbolic_points.py ifx
    python exe_manifolds.py ifx

### Figure 7

To generate the transmissivity data from Fig. 7(a), compile the `transmissivity_vs_N.f90` program as

    gfortran params_dp.f90 functions.f90 transmissivity_vs_N.f90 -fopenmp

or

    ifx params_dp.f90 functions.f90 transmissivity_vs_N.f90 -fopenmp

This program takes on five parameters: `num_ic`, `m`, `y0`, `esc_y`, and `N`. Run the program as

    ./a.out num_ic m y0 esc_y N

with `num_ic = 10000`, `m = -1`, and `N = 10000000000`. The `y0` and `esc_y` parameters have different values. For the blue curve, use `y0 = -1` and `esc_y = 1` and for the red curve, use `y0 = 1` and `esc_y = -1`.

To generate the rotation number data from Fig. 7(b), compile the `rotation_number.f90` program as

    gfortran params_dp.f90 functions.f90 rotation_number.f90 -fopenmp

or

    ifx params_dp.f90 functions.f90 rotation_number.f90 -fopenmp

and simply execute it as

    ./a.out

## Plotting the figures

After executing all the programs, simply run all cells in the Jupyter notebook `Plots.ipynb` to generate all the figures shown in the paper. The figures will be stored in a directory called `Figures/` which is created automatically by the Jupyter notebook.

## Citation

If you use this repository or parts of it in your work, please cite:

**Transport mechanisms associated with non-integer wavenumbers in a nontwist map**, *M. Rolim Sales et al.*

## Contact

For questions or feedback, feel free to [email me](mailto:rolim.sales.m@gmail.com).