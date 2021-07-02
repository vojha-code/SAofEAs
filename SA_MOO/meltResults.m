%% Insert this in the folder with the various results
% root=cd("Data");
contents = dir;
nEntries = 2;
folders=size(contents);
toSave = cell(10*folders(1)+1,9);

toSave(1,1) = {'Analysis'};
toSave(1,2) = {'Algorithm'};
toSave(1,3) = {'Function'};
toSave(1,4) = {'Metrics'};
toSave(1,5) = {'Parameter'};
toSave(1,6) = {'Mi/Si'};
toSave(1,7) = {'Sigma/Sti'};
toSave(1,8) = {'Mi/Si_Normalized'};
toSave(1,9) = {'Sigma/Sti_Normalized'};
Parameters(1).labels = {'lambda', 'mu-lambda_ratio', 'sigma0', 'alpha_mu','rescale_sigma0'};
Parameters(2).labels = {'pop_size', 'base_vector', 'diff_vector-pop_size_ratio', 'crossover_method', 'crossover_prob', 'beta_min','beta_ext'};
for i=1:size(contents)
    result = contents(i).name;
    if(startsWith(result,"20"))
        splitted = split(result,' ');
        splitted = split(splitted(3),'-');
        Analysis = splitted(1);
        Algorithm = splitted(2);
        Function = splitted(3);
        if (Algorithm== "CMA")
           Algorithm{1}='CMA-ES';
           Function = splitted(4);
        end
        folderContent = dir(strcat(result,"/ModelResults"));
        for j=1:size(folderContent)
            file = folderContent(j).name;
            if(endsWith(file,'.mat'))
                splitted = split(file,'.');
                Metrics=splitted(1);
                load(strcat(result,"/ModelResults/",file));
                if(Analysis=="SOBOL")
                    MiSi = Si;
                    SigmaSti = STi;
                else
                    MiSi = mi;
                    SigmaSti = sigma;
                end
                MiSi_Norm = (MiSi-min(MiSi))/(max(MiSi)-min(MiSi));
                SigmaSti_Norm = (SigmaSti-min(SigmaSti))/(max(SigmaSti)-min(SigmaSti));
                Size = size(MiSi);
                for p = 1:Size(2)
                    toSave(nEntries,1) = Analysis;
                    toSave(nEntries,2) = Algorithm;
                    toSave(nEntries,3) = Function;
                    toSave(nEntries,4) = Metrics;
                    if Algorithm=="CMA-ES"
                        toSave(nEntries,5) = Parameters(1).labels(p);
                    end
                    if Algorithm=="DE"
                        toSave(nEntries,5) = Parameters(2).labels(p);
                    end
                    toSave(nEntries,6) = {MiSi(p)};
                    toSave(nEntries,7) = {SigmaSti(p)};
                    toSave(nEntries,8) = {MiSi_Norm(p)};
                    toSave(nEntries,9) = {SigmaSti_Norm(p)};
                    nEntries = nEntries+1;
                end
            end
        end
    end
end
toSave=toSave(1:nEntries-1,:);
%cd(root);
cell2csv('Data/Results.csv',toSave)