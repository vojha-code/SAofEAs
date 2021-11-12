# Sensitivity Analysis Evolutionary Algorithms

The Sensitivity Analysis of Evolutinary Algorithms code reposiotry provides an comprehesive framework to study the influence of EAs hyper-paramters. This code reposetory builds on two senstivity analysis measures: elementry effect and varrience based effect.

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

We present a comprehensive global sensitivity analysis of two widely used single-objective and two multi-objective global optimization evolutionary algorithms (EAs). This work investigates the quality of influence hyperparameters have on the performances of EAs. Our methodology involves four EAs: covariance matrix adaptation evolutionary strategy (CMAES), differential evolution (DE), non-dominated sorting genetic algorithm III (NSGA-III), and multi-objective evolutionary algorithm based on decomposition (MOEA/D). We applied two sensitivity analysis methods, Morris and Sobol, to systematically analyze tunable hyperparameters of EAs on 23 single-objective functions and 10 multi-objective functions. Also, we studied the tendencies of hyperparameters’ direct effect and interaction effect with other hyperparameters. We present them as a comparative matrix where the diagonal from low direct effect and low interaction to high direct effect and high interaction shows the order and ranking of the hyperparameters. However, high interaction (or total effect) is a crucial measure of importance. Our investigation supports the research and innovation toward making different versions of DE as the type of DE algorithms emerged as one of the most influential hyperparameters for DE. It also supports the innovation in the type of strategy used for multi-objective decomposition in MOEA/D. Moreover, our results suggest the importance of correctly setting the initial step size for population generation in CMAES and crossover probability setting in NSGA-III.
