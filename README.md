# Semi-parametric Estimation and Simulation for Max-Mixture Models

This repository accompanies the methodology proposed in the paper titled: A semi-parametric estimation procedure for max-mixture spatial processes. It contains all the necessary R scripts and C files required to simulate, estimate, and validate various max-mixture spatial models.

---

## 📁 Structure

- `Estimation_and_Simulation_Workflow_commented.R`: Main workflow to simulate and estimate parameters
- `generate_TEG_process.R`: Contains the simulation function for generating data from the TEG process and Max-mixture process.
- `estimation_code.R`: Estimation procedures including custom likelihood and F-madogram-based estimators.
- `box_bar_density_plots.R`: Plotting utilities for RMSE boxplots, density plots, and bar plots.
- `*.c`, `*.o`, `*.so`, `f2c.h`: Supporting C code for efficient computation.

---

## 🔧 How to Use

### 1. **Load the project in RStudio**
Set working directory.

### 2. **Compile C source files**

From the R console:

```r
R CMD SHLIB src/cgj_new.c src/biv-nt.c
dyn.load("libs/cgj_new.so")
```

### 3. **Run the full estimation & simulation workflow**

```r
source("R/Estimation_and_Simulation_Workflow.R")
```

This script:
- Simulates spatial extreme data (e.g., TEG and max-mixture processes)
- Estimates parameters using semi-parametric and liklihood methods
- Computes RMSE and plots estimation performance

---

## 📈 Visualizations

`source(box_bar_density_plots.R)` generates:
- **Bar plots**: RMSE for each parameter  
- **Box plots**: Spread of estimated parameters  
- **Density plots**: Error distribution for each method

Each plot distinguishes between least squares and likelihood-based estimators.

---

## 📌 Methods

The estimation is based on:
- Least squares estimation via the proposed F-madogram 
- Simulation of spatial processes with known extremal structure

More details are in the paper (see citation below).

---

## 📜 License and Data Use

This code is released for academic, non-commercial use.

## 🧾 Citation

Please cite the following work if you use this code:

Ahmed, M. H., Maume-Deschamps, V., & Ribereau, P. (2025).  
*Semi-parametric estimation of extremal dependence structures with application to spatial extremes.*

---

## 👨‍🔬 Authors

- **Manaf H. Ahmed** – University of Mosul  
- **Véronique Maume-Deschamps** – Université Claude Bernard Lyon 1  
- **Pierre Ribereau** – Université Claude Bernard Lyon 1

📩 For questions: manaf.ahmed@uomosul.edu.iq
