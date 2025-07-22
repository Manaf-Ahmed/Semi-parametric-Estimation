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

### 3. Load the compiled shared objects  
```r
dyn.load("cgj_new.so")
dyn.load("rdavghol.so")  
```

### 4. Run the main workflow script  
```r
source("generate_TEG_process.R")
source("estimation_code.R")
source("Estimation_and_Simulation_Workflow.R")
```

These scripts:
- Simulate spatial extremes (e.g., TEG, max-mixture)
- Estimate parameters using F-madogram and censored likelihood
- All activities in the estimation, sImulation and plotting performance 

---

## 📈 Visualization

Running:
```r
source("box_bar_density_plots.R")
```
Generates:

- **Bar plots** – RMSE for each estimated parameter  
- **Box plots** – Variability of the estimated parameters  
- **Density plots** – Error distribution for each estimation method  

Visual comparisons are provided between **Least Squares (LS)** and **Likelihood** estimators.

---

## 📌 Methodology Summary

This repository implements:
- A novel **semi-parametric estimator based on the F-madogram**, compared with a parametric estimator using the **censored likelihood method** for estimating extremal dependence in spatial max-mixture processes.
- Simulated processes from max-stable and max-mixture models.  
- Estimation under different dependence structures: Smith (SM), Schlather (SC), and Brown–Resnick (BR).

Details are provided in the referenced paper below.

---

## 🌐 Data Availability

This code is provided for **academic and non-commercial use only**.  
Due to licensing restrictions from the Australian Bureau of Meteorology ([BoM](http://www.bom.gov.au)), original rainfall datasets are **not included**.  
Refer to [BoM data policy](http://www.bom.gov.au/other/copyright.shtml) for full usage terms.

**Station numbers used in the analysis:**

```
040020, 040051, 040094, 040138, 041025, 041128, 041310, 
053033, 055002, 055017, 055049, 060013, 061130, 062021, 
064009, 065000, 072056, 073007, 073110, 084015, 084025, 
084096, 084122, 035205, 044001, 040013, 043052, 052067, 
048027, 048025, 048176, 051049, 036040, 052003, 050108, 
060085, 058019, 066013, 062003
```

---

## 🧾 Citation

If you use this code, please cite the following work:

*Semi-parametric estimation of extremal dependence structures with application to spatial extremes.*

---

## 👨‍🔬 Authors

- **Manaf Ahmed** – Department of Statistics and Informatics, University of Mosul  
- **Véronique Maume-Deschamps** – ICJ, Université Claude Bernard Lyon 1  
- **Pierre Ribereau** – ICJ, Université Claude Bernard Lyon 1

📩 Contact: manaf.ahmed@uomosul.edu.iq
