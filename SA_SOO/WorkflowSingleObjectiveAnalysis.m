%{
Modified by V Ojha  - 21 April  2021

Earlier (2020) 
Modified by Alessio Greco from workflow_eet_hymod prepared by Francesca
Pianosi and Fanny Sarrazin
%} 
%% Step 1 (add paths)
my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))

%% Step 2a (setup Problem )
%global sd
%global ma;

load('Testbench.mat'); % problems_f variable is testbench
Problems = problems_f; % Problems
% Other data taken as gloabal variabls



%% Step 2b (setup Algorithms and Function to use)
Algorithms = {@() ypea_cmaes(), @() ypea_de()};
Metrics = {'BestSolution','ExecutionTime','BestSolutionGen'};
Termination = {'NEvaluations', 10000, 'MaxIter', 10000}; % MaxIter 10000  needed to avoid early stop for NEvaluations
myfun = 'evaluateSOEvAlgo' ; % Function that
parallel = true;

%% Step 2c (setup Parameters)
Parameters(1).labels = {'lambda', 'mu-lambda_ratio', 'sigma0', 'alpha_mu','rescale_sigma0'};
Parameters(1).Xmin = [10   , 0.1, 0.1, 0, 0];
Parameters(1).Xmax = [1000, 1  , 2  , 4, 1];
Parameters(2).labels = {'pop_size', 'base_vector', 'diff_vector-pop_size_ratio', 'crossover_method', 'crossover_prob', 'beta_min','beta_ext'};
Parameters(2).Xmin = [10 , 1, 0.01, 1, 0, 0, 0];
Parameters(2).Xmax = [1000, 5, 0.5 , 3, 1, 1, 2];
%% Step 3 (setup Analyses)

Analyses = {'Morris','MorrisLHS','SOBOL'};
design_types = {'trajectory', 'radial', ''}; % design_types for the analyses
NRun = 10;%0;

% Morris, MorrisLHS
r =  50;%r = 170; %

% Morris
L = 10;%L = 30  ; % number of levels in the uniform grid

% MorrisLHS, SOBOL
SampStrategy = 'lhs' ; % Latin Hypercube
% 'lhs': latin hypercube sampling
% 'rsu': random uniform
DistrFun = 'unif'; % Distribution function to use
%      'beta'  or 'Beta',
%      'bino'  or 'Binomial',
%      'chi2'  or 'Chisquare',
%      'exp'   or 'Exponential',
%      'ev'    or 'Extreme Value',
%      'f'     or 'F',
%      'gam'   or 'Gamma',
%      'gev'   or 'Generalized Extreme Value',
%      'gp'    or 'Generalized Pareto',
%      'geo'   or 'Geometric',
%      'hyge'  or 'Hypergeometric',
%      'logn'  or 'Lognormal',
%      'nbin'  or 'Negative Binomial',
%      'ncf'   or 'Noncentral F',
%      'nct'   or 'Noncentral t',
%      'ncx2'  or 'Noncentral Chi-square',
%      'norm'  or 'Normal',
%      'poiss' or 'Poisson',
%      'rayl'  or 'Rayleigh',
%      't'     or 'T',
%      'unif'  or 'Uniform',
%      'unid'  or 'Discrete Uniform',
%      'wbl'   or 'Weibull'

% SOBOL
SOBOL_weigth = 3; % How much times the VBSA should take more than the EET Analysis?

%% Step 4 (run the Various Analysis)
progressbar('buh','reset');
%progressbar('Analysis',numel(Analyses)-2);
for a=2:numel(Analyses)-1
    Analysis = Analyses(a);
    progressbar('Algotithm',numel(Algorithms));
    for n=2:numel(Algorithms)
        Algorithm = Algorithms{n};
        xmin=Parameters(n).Xmin;
        xmax=Parameters(n).Xmax;
        M=length(Parameters(n).labels);
        DistrPar = cell(M,1);
            
        fprintf('Algorithm %s \n', Algorithm().short_name); 
        for i=1:M
            DistrPar{i} = [ xmin(i) xmax(i)] ;
        end
        if(Analysis=="Morris")
            % use the sampling method originally proposed by Morris (1991):
            X = Morris_sampling(r,xmin,xmax,L); % (r*(M+1),M)
            if n == 1
                %save('X_Moris1.mat');
                MAT_M = load('X_Moris1.mat');
            else
                %save('X_Moris2.mat');
                MAT_M = load('X_Moris2.mat');
            end
            X = MAT_M.X;
            fprintf('Morris X:  %d %d \n', size(X)); 
        end
        if(Analysis=="MorrisLHS")
            % Latin Hypercube sampling strategy
            X = OAT_sampling(r,M,DistrFun,DistrPar,SampStrategy,design_types{a});
            if n == 1
                %save('X_LHS1.mat');
                MAT_LHS = load('X_LHS1.mat');
            else
                %save('X_LHS2.mat');
                MAT_LHS = load('X_LHS2.mat');
            end
            X = MAT_LHS.X;
            fprintf('MorrisLHS X:  %d %d \n', size(X)); 
        end
        
        if(Analysis=="SOBOL")
            % Sample parameter space using the resampling strategy proposed by
            % (Saltelli, 2008; for reference and more details, see help of functions
            % vbsa_resampling and vbsa_indices)
            N = 30;
            
            % Comment: the base sample size N is not the actual number of input
            % samples that will be evaluated. In fact, because of the resampling
            % strategy, the total number of model evaluations to compute the two
            % variance-based indices is equal to N*(M+2)
            X_AAT = AAT_sampling(SampStrategy,M,DistrFun,DistrPar,2*N);
            [ X, XB, XC ] = vbsa_resampling(X_AAT) ;
            if n == 1
                %save('XA_SOBOL1.mat');
                %save('XB_SOBOL1.mat');
                %save('XC_SOBOL1.mat');

                MAT_A = load('XA_SOBOL1.mat');
                MAT_B = load('XB_SOBOL1.mat');
                MAT_C = load('XC_SOBOL1.mat');
            else
                %save('XA_SOBOL2.mat');
                %save('XB_SOBOL2.mat');
                %save('XC_SOBOL2.mat');

                MAT_A = load('XA_SOBOL2.mat');
                MAT_B = load('XB_SOBOL2.mat');
                MAT_C = load('XC_SOBOL2.mat');
            end
            

            X = MAT_A.X;
            XB = MAT_B.XB;
            XC = MAT_C.XC;
            
            fprintf('SOBOL X:  %d %d \n', size(X)); 
            fprintf('SOBOL XB: %d %d \n', size(XB)); 
            fprintf('SOBOL XC: %d %d \n', size(XC)); 
        end
        %break
        %error('stop');
        
        progressbar('Problem',numel(Problems));
        for i=1:numel(Problems)
            if(i == 8)% && (i ~= 19) && (i ~= 20)
                fprintf('Runing problem...: %d \n',i);
                
                Problem=Problems(i);
                AlgoObj = Algorithm();
                folder = strcat(datestr(now, 'yyyy-mm-dd HH-MM-SS')," ",Analysis,"-",AlgoObj.short_name,"-Problem",int2str(i));
                
                if(exist(strcat(folder,'/ModelEval'),'dir')~=7)
                    mkdir(strcat(folder,'/ModelEval'));
                end
                root = cd(strcat(folder,'/ModelEval'));
                %Execution
                Y = model_evaluation(myfun,X, Algorithm, Problem, Metrics, Termination, NRun, parallel) ; % size (r*(M+1),3)
                if (Analysis=="SOBOL")
                    YB = model_evaluation(myfun,XB, Algorithm, Problem, Metrics, Termination, NRun, parallel); % size (r*(M+1),3)
                    YC = model_evaluation(myfun,XC, Algorithm, Problem, Metrics, Termination, NRun, parallel); % size (r*(M+1),3)
                end
                cd(root);
                
                if(exist(strcat(folder,'/ModelResults'),'dir')~=7)
                    mkdir(strcat(folder,'/ModelResults'));
                end
                root = cd(strcat(folder,'/ModelResults'));
                
                for j=1:length(Y(1,:))
                    filename = Metrics{j};
                    if(Analysis=="SOBOL")
                        % Option 3: SOBOL
                        [ Si, STi ] = vbsa_indices(Y(:,j),YB(:,j),YC(:,j));
                        save(filename,'Si','STi');
                    else
                        % Option 1/2: EET
                        [ mi, sigma ] = EET_indices(r,xmin,xmax,X,Y(:,j),design_types{a});
                        save(filename,'mi','sigma');
                    end
                end
                cd(root);
                if(Analysis=="SOBOL")
                    save(['Data/',char(Analysis),'-',AlgoObj.short_name,'-Problem',int2str(i),'.mat'],'Y','YB','YC')
                else
                    save(['Data/',char(Analysis),'-',AlgoObj.short_name,'-Problem',int2str(i),'.mat'],'Y')
                end
                progressbar('Problem',1);
            else
                %disp('do not run ');
                fprintf('do not run problem: %d \n',i);
            end
        end
        progressbar('Algotithm',1);
    end
    %progressbar('Analysis',1);
end

