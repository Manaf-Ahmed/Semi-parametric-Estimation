# Semi-parametric Estimation and Simulation for Max-Mixture Models

This repository supports the methodology proposed in the paper titled:  
**"A Semi-parametric Estimation Procedure for Max-Mixture Spatial Processes"**.  
It includes all necessary R scripts and C source files to simulate, estimate, and validate spatial max-mixture models.

---

## 📁 Repository Structure

- `Estimation_and_Simulation_Workflow.R` – Main script for simulating and estimating model parameters  
- `generate_TEG_process.R` – Functions for simulating TEG and max-mixture processes  
- `estimation_code.R` – Semi-parametric and censored likelihood-based estimation functions  
- `box_bar_density_plots.R` – Functions to visualize RMSE, boxplots, and density plots  
- `*.c`, `f2c.h`, `*.so` – C code for fast numerical computations  

---

## 🔧 How to Use

### 1. Set your working directory in R
Make sure all scripts and compiled libraries are in the same folder.

### 2. Compile the C source files (from R console)
```r
system("R CMD SHLIB cgj_new.c biv-nt.c")
system("R CMD SHLIB rdavghol.c")  # If needed
```
### 3. Load the C functions  
```r
dyn.load("cgj_new.so")
dyn.load(paste("rdavghol.so")
```
### 4. Run the main workflow
```r
source("generate_TEG_process.R")
source("estimation_code.R")
source("Estimation_and_Simulation_Workflow.R")
```

These scripts:
- Simulate spatial extremes (e.g., TEG, max-mixture)
- Estimate parameters using F-madogram and likelihood
- Compute RMSE and generate performance plots

---

## 📈 Visualization

Running 
```r
source("box_bar_density_plots.R")
```
provides:

- **Bar plots** – RMSE for each estimated parameter  
- **Box plots** – Variability of the estimated parameters  
- **Density plots** – Error distribution for each estimation method  

Visual comparisons are shown between **Least Squares (LS)** and **Likelihood** estimators.

---

## 📌 Methodology Summary

This repository implements:
- A novel **semi-parametric estimator based on the F-madogram** is proposed and compared with a parametric estimator using the **censored likelihood method** for estimating extremal dependence in spatial max-mixture processes
- Simulated processes from max-stable and max-mixture models  
- Estimation under different dependence structures: Smith, Schlather, Brown–Resnick  

Detailed methodology can be found in the referenced paper below.

---

## 📜 License and Data Policy

This code is provided for **academic and non-commercial use only**.  
Due to licensing restrictions from the Australian Bureau of Meteorology (http://www.bom.gov.au), original rainfall datasets are **not included**.  
Refer to [BOM data policy](http://www.bom.gov.au/other/copyright.shtml) for usage details.

---

## 🧾 Citation

If you use this code, please cite the following work:

> Ahmed, M. H., Maume-Deschamps, V., & Ribereau, P. (2025).  
> *Semi-parametric estimation of extremal dependence structures with application to spatial extremes.*

---

## 👨‍🔬 Authors

- **Manaf  Ahmed** – Department of Statistics and Informatics, University of Mosul  
- **Véronique Maume-Deschamps** – ICJ, Université Claude Bernard Lyon 1  
- **Pierre Ribereau** – ICJ, Université Claude Bernard Lyon 1

📩 Contact: manaf.ahmed@uomosul.edu.iq
