    function Score = evaluateMOEvAlgo(Parameters, Algorithm, Problem, Metrics, Termination, NRun, parallel, FilePath)
    nMetrics=length(Metrics);
    result = zeros(NRun,nMetrics);
    fileName = strjoin(strsplit(num2str(Parameters)),'-');
    %mkdir(folder);
    %fprintf('Progress sample : %s \n',fileName);
    
    %fprintf('Progress:\n');
    %fprintf(['\n' repmat('.',1,NRun) '\n\n']);
    
    if(parallel)
        parfor i=1:NRun
            result(i,:) = run(Parameters, Algorithm, Problem, Metrics, Termination);
            %saveypea([folder,'/Run',int2str(i),'.mat'],ypea_algo)
            %fprintf(" |-> (%d)",i);
            %fprintf('\b|\n');
        end
    else
        for i=1:NRun
            result(i,:) = run(Parameters, Algorithm, Problem, Metrics, Termination);
            %save([folder,'/Run',int2str(i),'.mat'],'ypea_algo')
            %fprintf(" |-> (%d)",i);
            %error('stop');
        end
    end
    %cd ..
    %cd(strcat(FilePath,'/ModelEval'))
    
    writetable(array2table(result,'VariableNames',Metrics), [strcat(FilePath,fileName,"-metrics.csv")]);
    %result
    Score = mean(result,1);
    %error('stop');
    end

    function saveypea(name, var)
    save(name,'var');
    end

    function result=run(Parameters, Algorithm, Problem, Metrics, Termination)
    result=zeros(numel(Metrics),1);
    warning('off','all');
    
    emoResults =  [];
    %Parameters, Algorithm, Problem, Metrics, Termination
    if (func2str(Algorithm) == "MOEAD")
        %pop = round(Parameters(1)); % 100; % Population size  % numeric
        %pSBX = Parameters(2); %0.5;      % 1 crossover
        %sbxDI = round(Parameters(3)); %100;     % 2 crossover % numeric
        %pm = Parameters(4); %0.6;        % 3 mutation
        %pmDI = round(Parameters(5)); %150;      % 4 mutation  % numeric
        
        %mode = round(Parameters(6)); %1;        % 7 mode 
        %neighbour = Parameters(7); %0.3; % 8  neighbourhood   
        
        % We neeed dummay parameter
        tournament = 2;  % DUMMY
        emoResults = platemo('algorithm',{Algorithm, Parameters(2), round(Parameters(3)), Parameters(4), round(Parameters(5)), tournament, round(Parameters(6)), Parameters(7)},'problem',Problem{1,1},'M',Problem{1,2},'D',(Problem{1,2}+10-1),'N',round(Parameters(1)),'maxFE',Termination);
        %results = platemo('algorithm',{Algorithm, pSBX, sbxDI, pm, pmDI, tournament, mode, neighbour},'problem',Problem{1,1},'M',Problem{1,2},'D',(Problem{1,2}+10-1),'N',pop,'maxFE',Termination);
    else
        %NSGAIII 
        %pop = round(Parameters(1)); % 100; % Population size   % numeric
        %pSBX = Parameters(2); %0.5;      % 1 crossover
        %sbxDI = round(Parameters(3)); %100;     % 2 crossover  % numeric
        %pm = Parameters(4); %0.6;        % 3 mutation
        %pmDI = round(Parameters(5)); %150;      % 4 mutation   % numeric
        %tournament = round(Parameters(6)); % 5;  % 6 turnament % numeric

        % We neeed dummay parameter
        mode = 1;        % DUMMY
        neighbour = 0.3; % DUMMY
        emoResults = platemo('algorithm',{Algorithm, Parameters(2), round(Parameters(3)), Parameters(4), round(Parameters(5)), round(Parameters(6)), mode, neighbour},'problem',Problem{1,1},'M',Problem{1,2},'D',(Problem{1,2}+10-1),'N',round(Parameters(1)),'maxFE',Termination);
        %results = platemo('algorithm',{Algorithm, pSBX, sbxDI, pm, pmDI, tournament, mode, neighbour},'problem',Problem{1,1},'M',Problem{1,2},'D',(Problem{1,2}+10-1),'N',pop,'maxFE',Termination);
    end
    % emoResults = platemo('algorithm',{Algorithm, pSBX, sbxDI, pm, pmDI, tournament, mode, neighbour},'problem',Problem{1,1},'M',Problem{1,2},'D',(Problem{1,2}+10-1),'N',pop,'maxFE',Termination);

    gdITR = 0;
    igdITR = 0;
    hvITR = 0;
    for index = 1:length(emoResults)
        cur_result = emoResults{index,2};
        pro = Problem{1,1}('M',Problem{1,2},'D',(Problem{1,2}+10-1));
        gdITR  = gdITR + GD(cur_result,pro.optimum);
        igdITR  = igdITR + IGD(cur_result,pro.optimum);
        hvITR  = hvITR + HV(cur_result,pro.optimum); % Original HV implementation
        %hvITR  = hvITR + HV1(cur_result); % Custom HV implementation
        %fprintf("gd %d, idg %d, hv %d \n",gd,igd,hv);
    end
    gd = (gdITR/length(emoResults));
    igd = (igdITR/length(emoResults));
    hv = (hvITR/length(emoResults));
    %hv = HV(emoResults{end},pro.optimum);
    %fprintf("gd %d, idg %d, hv %d \n",gd,igd,hv);
    
    for j = 1:numel(Metrics)
        switch Metrics{j}
            case "GD" %Generational Distance
                result(j) = gd;
            case "IGD" %Inverse Generational Distance
                result(j) = igd;
            case "HV" %Hypervolume Indicator
                result(j) = hv;
        end
    end
    %result

    end