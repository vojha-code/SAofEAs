%{
V Ojha data prosessing for Moriss and SOBAL results
%}
%% Root dir
my_dir = pwd;
root = cd(my_dir);
addpath(genpath(my_dir));
fprintf('current directory %s',my_dir);

%% Check folders
folder = 'SA_MO_Results_Collection';
folder = fullfile(root, folder);
fprintf('\n directory investigated %s',folder);
%% Algorithm Param and Settings and Resulst collection
listing = dir(folder);

MOEAD = ["pop" , "pSBX", "sbxDI", "pm", "pmDI", "mode", "neighbour"];
NSGAIII = ["pop" , "pSBX" , "sbxDI" , "pm" , "pmDI" , "tournament"];

Algo = {'MOEAD', 'NSGAIII'};
Metric = ["GD";"HV";"IGD"];
for metric =1:length(Metric)
    for k=1:numel(Algo)
        algo= Algo(k)
        param = MOEAD;
        if string(algo) == "NSGAIII"
            param = NSGAIII;
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
                    modelResulstBestSol = load(fullfile(fullFolderModelResult,listFiles(2+metric).name));
                    fprintf('Model Results in %s accessed %s \n', nextFolder, listFiles(2+metric).name);
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
        writematrix(ResultsMoris,"ResultsMoris_"+string(algo)+"_"+Metric(metric)+"_20var.csv");
        writematrix(ResultsMorisLHS,"ResultsMorisLHS_"+string(algo)+"_"+Metric(metric)+"_20var.csv");
        writematrix(ResultsSOBAL,"ResultsSOBAL_"+string(algo)+"_"+Metric(metric)+"_28var.csv");
    end
end