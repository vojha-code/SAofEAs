function Score = evaluateSOEvAlgo(Parameters, Algorithm, Problem, Metrics, Termination, NRun, parallel)
    nMetrics=length(Metrics);
    result = zeros(NRun,nMetrics);
    folder = strjoin(strsplit(num2str(Parameters)),'-');
    %mkdir(folder);
    % Initialize algorithm
    ypea_algo = Algorithm();
    % Set problem
    ypea_problem = Problem;
    % Set parameters
    ypea_algo.SetParameters(Parameters);
    % Set termination
    ypea_algo.SetTermination(Termination);
    % Disable display
    ypea_algo.display=false;
    %fprintf('Progress sample : %s \n',folder);
    %fprintf(['\n' repmat('.',1,NRun) '\n\n']);
    if(parallel)
         parfor i=1:NRun
            result(i,:) = run(ypea_algo, ypea_problem, Metrics);
            %saveypea([folder,'/Run',int2str(i),'.mat'],ypea_algo)
            %fprintf('\b|\n');

        end       
    else
        for i=1:NRun
            result(i,:) = run(ypea_algo, ypea_problem, Metrics);
            %save([folder,'/Run',int2str(i),'.mat'],'ypea_algo')
            %fprintf('\b|\n');
            %error('quit');
            %quit(1);
        end
    end
    writetable(array2table(result,'VariableNames',Metrics), [strcat(folder,"-metrics.csv")]);
    Score = mean(result,1);
end

function saveypea(name, var)
    save(name,'var');
end

function result=run(ypea_algo, ypea_problem, Metrics)
    result=zeros(numel(Metrics),1);
    warning('off','all');
    ypea_algo.solve(ypea_problem);
    for j = 1:numel(Metrics)
        switch Metrics{j}
            case "BestSolution"
                result(j) = ypea_algo.best_obj_value;
            case "ExecutionTime"
                result(j) = ypea_algo.run_time;
            case "BestSolutionGen"
                best_sol_gen =  find(ypea_algo.best_obj_value_history == ypea_algo.best_obj_value,1);
                if(isempty(best_sol_gen)) 
                    result(j) = 1;
                else
                    result(j) = best_sol_gen;
                end
        end
    end
end