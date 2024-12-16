# Long-term persistence of transposable elements activity
This repository contains the code and data files required to reproduce the figures and simulation results presented in the paper "Maintenance of long-term transposable element activity through regulation by nonautonomous elements" by Omole, A. D., and Czuppon, P.

## Reproducing Figures from the Manuscript

The following instructions detail how to reproduce Figures 2, 3, 4, and 5 from the manuscript using the provided code and scripts:

### Figure 2
- **Folder**: `TEs_stochastic_simulation`
- **Data Generation**:  
  Run the Julia script `pop_TE_dynamics.jl` to generate the dataset `pop_TE_dynamics.csv`.
- **Data Processing**:  
  Process the dataset using the Python script `Pop_TE_dynamics.py` to create the figure.

### Figures 3 and 5
- **Folder**: `Stationary_distribution`
- **Data Generation**:  
  Run the Julia script `Pop_stationary_distribution.jl` to generate the required dataset.
- **Data Processing**:  
  Use the Python script `Pop_stationary_distribution.py` to process the generated data and produce the figures.

### Figure 4
- **Folder**: `Relative_mean_and_variance`
- **Data Generation**:  
  Run the Julia script `Data_variance_alphavalues_proposition2.jl` to generate the necessary dataset.
- **Data Processing**:  
  Process the data using the Python script `Variance_Alphavalues_proposition2.py` to create the figure.

---

Each folder contains the relevant scripts to reproduce the figures. Ensure all dependencies for Julia and Python are installed before running the scripts.
