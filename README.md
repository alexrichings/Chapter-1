# Chapter-1

This repository contains the scripts supporting the work titled: 

**"Using finite mixture models to leverage multiplex immunoassay data from hyperendemic settings"**

This study investigates the use of boosts in antibody levels as a metric of exposure instead of seroconversion to combat identifiability issues in hyperendemic settings, and is applied to large-scale two-sample longitudinal serological data collected between 2021 and 2022 in the Dominican Republic. 

----------------------------------------------------------------------------------------------------------------

# Repository Structure 

----------------------------------------------------------------------------------------------------------------

# Script Overview 

## Dominican Republic (DR) scripts 

----------------------------------------------------------------------------------------------------------------

### 01-setup_inspect
*Purpose:*
  - Import raw serological data
  - Clean variables
  - Visual inspection of data 

*Output:*
  - Paired serological data file (df_wide)

----------------------------------------------------------------------------------------------------------------

### 02-fmm-dr
*Purpose:*
  - Run two-component finite mixture model (FMM) in a frequentist framework 

*Output:*
  - Indivdual probability of boosting for each person per dengue serotype (inputs_dr)




  - FMM parameters
  - FMM component plots:
    - deng1_fmm_prob_plot.png
    - deng2_fmm_prob_plot.png
    - deng3_fmm_prob_plot.png
    - deng4_fmm_prob_plot.png
    - facet_fmm_plot.png
   
----------------------------------------------------------------------------------------------------------------

### 03-boost-regression

