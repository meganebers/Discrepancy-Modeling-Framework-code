Matlab code to accompany the paper submitted to SIADS titled "Discrepancy Modeling Framework: Learning missing physics, modeling systematic residuals, and disambiguating between deterministic and random effects" by Megan R. Ebers, Katherine M. Steele, and J. Nathan Kutz (University of Washington 2022). The preprint can be found at https://arxiv.org/abs/2203.05164.

Generate_All_Data generates, trains, and tests 4 data-driven models (GPR, opt-DMD, SINDy, and feedforward NN) across 3 levels of additive Gaussian noise for both types of discrepanices: missing physics/dynamics and systematic state-space residuals. We evaluate this framework for the Van der Pol oscillator, the Lorenz 63 attractor, and the Burgers wave equation.

Process_Data_Results generates figures based on the results from Generate_All_Data. 

Please direct all questions to mebers@uw.edu