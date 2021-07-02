%{
V Ojha data prosessing for Moriss and SOBAL results
%}
%% Root dir
my_dir = pwd;
root = cd(my_dir);
addpath(genpath(my_dir));
fprintf('current directory %s',my_dir);

%% Check folders
folder = 'SO_Results_NF_10K';
folder = fullfile(root, folder);
fprintf('directory investigated %s',folder);
%% Algorithm Param and Settings and Resulst collection
listing = dir(folder);

CMAES = ["lambda", "mu-lambda_ratio", "sigma0", "alpha_mu", "rescale_sigma0"];
DE = ["pop_size", "base_vector", "diff_vector-pop_size_ratio", "crossover_method", "crossover_prob", "beta_min", "beta_ext"];

Algo = {'DE', 'CMA-ES'};
for k=1:numel(Algo)
    algo= Algo(k)
    param = CMAES;
    if string(algo) == "DE"
        param = DE;    
    end

    ResultsMoris = [];
    ResultsMorisLHS  = [];
    ResultsSOBAL = [];
    ResultsMoris = [ResultsMoris; param, param];
    ResultsMorisLHS  = [ResultsMorisLHS; param, param];
    ResultsSOBAL = [ResultsSOBAL; param, param];
    for i = 1:length(listing)
        nextFolder = listing(i).name;
        fullFolderPath = fullfile(folder, nextFolder);
        if strlength(nextFolder) > 5
            fullFolderModelResult = fullfile(fullFolderPath, 'ModelResults');
            listFiles = dir(fullFolderModelResult);
            if contains(nextFolder, string(algo))
                modelResulstBestSol = load(fullfile(fullFolderModelResult,listFiles(3).name));
                fprintf('Model Results in %s accessed %s \n', nextFolder, listFiles(3).name);
                if contains(nextFolder,'Morris-')
                    ResultsMoris =  [ResultsMoris; modelResulstBestSol.mi, modelResulstBestSol.sigma];
                end
                if contains(nextFolder,'MorrisLHS-')
                    ResultsMorisLHS =  [ResultsMorisLHS; modelResulstBestSol.mi, modelResulstBestSol.sigma];
                end
                if contains(nextFolder,'SOBOL-')
                    ResultsSOBAL =  [ResultsSOBAL; modelResulstBestSol.STi, modelResulstBestSol.Si];
                end
            end
        end
        % all of your actions for filtering and plotting go here
    end
    writematrix(ResultsMoris,'ResultsMoris_'+string(algo)+'_Sample1K_Term10K.csv');
    writematrix(ResultsMorisLHS,'ResultsMorisLHS_'+string(algo)+'_Sample1K_Term10K.csv');
    writematrix(ResultsSOBAL,'ResultsSOBAL_'+string(algo)+'_Sample1K_Term10K.csv');
end