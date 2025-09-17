# Interval-Censored Data Analysis under the Cox Proportional Hazards Model

This project provides a comprehensive study of **interval-censored data analysis under the Cox proportional hazards model**.

## Functions

The functions in **`functions.r`** support the following tasks:

- Estimating the **sieve spline-based maximum likelihood estimator (MLE)** of the Cox regression parameters and the hazard function.  
- Estimating the variance of the sieve MLE of the Cox regression parameters using two established methods.  
- Generating knot sequences for **B-splines or I-splines**.  

## Scripts

### `simulation.r`
Compares **our proposed likelihood ratio testing approach** with existing **Wald testing approaches**.  

### `hemophilia.r`
Applies different testing approaches to analyze a publicly available dataset.  

### `UMvsEM.r`
Compares computing time and estimation bias between **the proposed unconstrained maximization (UM)** method and **the EM algorithm implemented in the R package ICsurv**.
