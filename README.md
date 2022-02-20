# Sensitivity Analysis Evolutionary Algorithms

The Sensitivity Analysis of Evolutinary Algorithms code reposiotry provides an comprehesive framework to study the influence of EAs hyper-paramters. This code reposetory builds on two senstivity analysis measures: elementry effect and varrience based effect.

Prepreint Submitted to Swarm and Evolutionary computation, Elsevier

All code scriipts are in MATLAB. The are optimized for MATLAB version 2020. 
The respective libraries used should be refered to to theire original author/version.  

## Evolutionary Algorithms Studied
- Single Objective Algorithms (SA_SOO folder)
    - Diffrential Evolution (DE)
    - Covariance Matrix Adaptation Evolution Strategy (CMA-ES)
- Multi-Ojective Algorithms (SA_MOO folder)
    - Non-dominated Sorting Genetic Algorithm III (NSGA-III)
    - Multiobjective Evolutionary Algorithm Based on Decomposition (MOEA/D)
- Results of Sansitivity Analysis of Evolutionary Algorithms (SA_EA_Results)
    - Results MOO
    - Results SOO
    - Py Processess Data (code for processing data)
- README


## Results

We present a comprehensive global sensitivity analysis of two single-objective and two multi-objective state-of-the-art global optimization evolutionary algorithms as an algorithm configuration problem . That is, we investigate the quality of influence hyperparameters have on the performance of algorithms in terms of their direct effect
and interaction effect with other hyperparameters. Using three sensitivity analysis methods, Morris LHS, Morris, and Sobol, to systematically analyze tunable hyperparameters of covariance matrix adaptation evolutionary strategy, differential evolution, non-dominated sorting genetic algorithm III, and multi-objective evolutionary algorithm based on decomposition, the framework reveals the behaviors of hyperparameters to sampling methods and performance metrics. That is, it answers questions like what hyperparameters influence patterns, how they interact, how much they interact, and how much their direct influence is. Consequently, the ranking of hyperparameters suggests their order of tuning, and the pattern of influence reveals the stability of the algorithms.

<img src="https://github.com/vojha-code/SAofEAs/blob/master/SA_EA_Results/DE_param.png" alt="DE and CMA-ES" width="700">
<img src="https://github.com/vojha-code/SAofEAs/blob/master/SA_EA_Results/NSGA_III.png" alt="DE and CMA-ES" width="700">
