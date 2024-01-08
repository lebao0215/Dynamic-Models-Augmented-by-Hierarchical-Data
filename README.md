# Dynamic-Models-Augmented-by-Hierarchical-Data
Dynamic Models Augmented by Hierarchical Data: An Application Of Estimating HIV Epidemics At Sub-National And Sub-Population Level
ThaiMain.R is the main code.
First, set the working directory to the same location as ThaiMain.R, and then install the following packages:
"brms", "bayesplot", "splines"
KL.R is called by ThaiMain.R. It uses KL divergence to select the optimal sample size for the augmented data.
ThailandEPP.R runs the EPP model. We need to specify the site from 1 to 4, and specify the input data with or without the augmented data.
