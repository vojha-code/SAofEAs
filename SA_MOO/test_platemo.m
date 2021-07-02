my_dir = pwd ; % use the 'pwd' command if you have already setup the Matlab
% current directory to the SAFE directory. Otherwise, you may define
% 'my_dir' manually by giving the path to the SAFE directory, e.g.:
% my_dir = '/Users/francescapianosi/Documents/safe_R1.0';

% Set current directory to 'my_dir' and add path to sub-folders:
cd(my_dir)
addpath(genpath(my_dir))
%%

warning('off');

Problem = { @ZDT1,
    @ZDT2,
    @ZDT3,
    @ZDT4,
    @ZDT5,
    @ZDT6,
    @DTLZ1,
    @DTLZ2,
    @DTLZ3,
    @DTLZ4,
    @DTLZ5,
    @DTLZ6,
    @DTLZ7,
    @DTLZ8,
    @DTLZ9,
    @C1_DTLZ1,
    @C1_DTLZ3,
    @C2_DTLZ2,
    @C3_DTLZ4,
    @CDTLZ2,
    @DC1_DTLZ1,
    @DC1_DTLZ3,
    @DC2_DTLZ1,
    @DC2_DTLZ3,
    @DC3_DTLZ1,
    @DC3_DTLZ3,
    @IDTLZ1,
    @IDTLZ2,
    @SDTLZ1,
    @SDTLZ2,
    @WFG1,
    @WFG2,
    @WFG3,
    @WFG4,
    @WFG5,
    @WFG6,
    @WFG7,
    @WFG8,
    @WFG9
    };
%%
Problems = [ {@DTLZ1,3}; {@DTLZ1,5};  {@DTLZ1,8};  {@DTLZ1,10}; {@DTLZ1,15}; % K. Deb, L. Thiele, M. Laumanns and E. Zitzler (DTLZ)
    %{@DTLZ2,3}; {@DTLZ2,5};  {@DTLZ2,8};  {@DTLZ2,10}; {@DTLZ2,15};
    %{@DTLZ3,3}; {@DTLZ3,5};  {@DTLZ3,8};  {@DTLZ3,10}; {@DTLZ3,15};
    %{@DTLZ4,3}; {@DTLZ4,5};  {@DTLZ4,8};  {@DTLZ4,10}; {@DTLZ4,15};
    %{@DTLZ5,3} {@DTLZ5,5}  {@DTLZ5,8}  {@DTLZ5,10} {@DTLZ5,15}
    %{@DTLZ6,3} {@DTLZ6,5}  {@DTLZ6,8}  {@DTLZ6,10} {@DTLZ6,15}
    %{@DTLZ7,3} {@DTLZ7,5}  {@DTLZ7,8}  {@DTLZ7,10} {@DTLZ7,15}
    %{@DTLZ8,3} {@DTLZ8,5}  {@DTLZ8,8}  {@DTLZ8,10} {@DTLZ8,15}
    %{@DTLZ9,3} {@DTLZ9,5}  {@DTLZ9,8}  {@DTLZ9,10} {@DTLZ9,15}
    {@CDTLZ2,3}; {@CDTLZ2,5};  {@CDTLZ2,8};  {@CDTLZ2,10}; {@CDTLZ2,15}; % Convex DTLZ2
    %{@C1_DTLZ1,3}; {@C1_DTLZ1,5};  {@C1_DTLZ1,8};  {@C1_DTLZ1,10}; {@C1_DTLZ1,15}; %  Constrained DTLZ1 % Jain, H. and K. Deb. "An Evolutionary Many-Objective Optimization Algorithm Using Reference-Point-Based Nondominated Sorting Approach, Part II: Handling Constraints and Extending to an Adaptive Approach." EEE Transactions on Evolutionary Computation, 18(4):602-622, 2014
    %{@C1_DTLZ3,3} {@C1_DTLZ3,5}  {@C1_DTLZ3,8}  {@C1_DTLZ3,10} {@C1_DTLZ3,15},
    %{@C2_DTLZ2,3} {@C2_DTLZ2,5}  {@C2_DTLZ2,8}  {@C2_DTLZ2,10} {@C2_DTLZ2,15},
    %{@C3_DTLZ4,3} {@C3_DTLZ4,5}  {@C3_DTLZ4,8}  {@C3_DTLZ4,10} {@C3_DTLZ4,15},
    {@IDTLZ1,3}; {@IDTLZ1,5};  {@IDTLZ1,8}; {@IDTLZ1,10}; {@IDTLZ1,15};   % Inverted DTLZ1
    {@IDTLZ2,3}; {@IDTLZ2,5};  {@IDTLZ2,8}; {@IDTLZ2,10}; {@IDTLZ2,15};   % Inverted DTLZ2
    %{@WFG1,3} {@WFG1,5}  {@WFG1,8}  {@WFG1,10} {@WFG1,15}
    %{@WFG2,3} {@WFG2,5}  {@WFG2,8}  {@WFG2,10} {@WFG2,15}
    {@WFG3,3}; {@WFG3,5};  {@WFG3,8};  {@WFG3,10}; {@WFG3,15};
    %{@WFG4,3} {@WFG4,5}  {@WFG4,8}  {@WFG4,10} {@WFG4,15}
    %{@WFG5,3} {@WFG5,5}  {@WFG5,8}  {@WFG5,10} {@WFG5,15}
    {@WFG6,3}; {@WFG6,5};  {@WFG6,8};  {@WFG6,10}; {@WFG6,15};
    {@WFG7,3}; {@WFG7,5};  {@WFG7,8};  {@WFG7,10}; {@WFG7,15};
    %{@WFG8,3} {@WFG8,5}  {@WFG8,8}  {@WFG8,10} {@WFG8,15}
    %{@WFG9,3} {@WFG9,5}  {@WFG9,8}  {@WFG9,10} {@WFG9,15}
    ];

%%
Algorithms = {@MOEAD, @NSGAIII};

Objectives = [3, 5, 8, 10, 15];  % Number of objectives
Dimendsion = 20; % Number of variables
save = 5;
maxFE = 10000;

algo = 2;
problem = 1;
obj = 5;

Parameters =[ 20    ,   0.1,   2,     0.1,   2,     1,     0.1];
pop = round(Parameters(1)); % 100; % Population size  % numeric
pSBX = Parameters(2); %0.5;      % 1 crossover
sbxDI = round(Parameters(3)); %100;     % 2 crossover % numeric
pm = Parameters(4); %0.6;        % 3 mutation
pmDI = round(Parameters(5)); %150;      % 4 mutation  % numeric

mode = round(Parameters(6)); %1;        % 7 mode 
neighbour = Parameters(7); %0.3; % 8  neighbourhood   

% We neeed dummay parameter
tournament = 2;  % DUMMY
%%
igd = [];
gd = [];
hv = [];

%[Dec,Obj,Con] = platemo('algorithm',@MOEAD,'problem',@ZDT1,'N',N,'M',M,'D',D,'save',save,'maxFE',maxFE);

for algo = 1:length(Algorithms)
    for problem =1:length(Problems)
        
        fprintf("%s, %s (%d)\n",func2str(Algorithms{algo}), func2str(Problems{problem,1}),Problems{problem,2});
        
        %results = platemo('algorithm',Algorithm{2},'problem',Problem{1},'M',5,'N',N,'D',D,'save',save,'maxFE',maxFE);
        results = platemo('algorithm',{Algorithms{algo}, pSBX, sbxDI, pm, pmDI, tournament, mode, neighbour},'problem',Problems{problem,1},'N',pop,'M',Problems{problem,2},'D',(Problems{problem,2}+10-1),'maxFE',maxFE);
        %results = platemo('algorithm',{Algorithm{algo}, pSBX, sbxDI, pm, pmDI, offsping, tournament, mode, neighbour},'problem',Problem{problem,1},'M',Problem{problem,2},'D',(Problem{problem,2}+10-1),'N',N,'maxFE',maxFE);
        
        %platemo('algorithm',@NSGAII,'problem',@DTLZ2,'M',5,'D',40,'maxFE',20000,'save',10);
        
        %platemo('algorithm',@MOEAD,'problem',@ZDT1,'N',100,'M',2,'D',10,'save',10,'maxFE',(100*1000))
        
        
        %fprintf(" \n Outputs: ");
        gdITR = 0;
        igdITR = 0;
        hvITR = 0;
        for index = 1:length(results)
            cur_result = results{index,2}; % cur_result = results{end};
            pro = Problems{problem,1}('M',Problems{problem,2},'D',(Problems{problem,2}+10-1));
            gdITR  = gdITR + GD(cur_result,pro.optimum);
            igdITR  = igdITR + IGD(cur_result,pro.optimum);
            %hvITR  = hvITR + HV(cur_result,pro.optimum);
            hvITR  = hvITR + HV1(cur_result); % Custome HV 
            %fprintf("gd %d, idg %d, hv %d \n",gd,igd,hv);
        end
        
        gd  = [gd; (gdITR/length(results))];
        igd  = [igd; (igdITR/length(results))];
        hv  = [hv; (hvITR/length(results))];
        %fprintf(" gd %d, idg %d, hv %d \n",gd,igd,hv);
        %}
        
    end
end

%%
fprintf("\nResults: \n");
i = 1;
for algo = 1:length(Algorithms)
    for problem =1:length(Problems)
        fprintf("%s, %s (%d):  gd %d, idg %d, hv %d \n",func2str(Algorithms{algo}), func2str(Problems{problem,1}), Problems{problem,2}, gd(i),igd(i),hv(i));
        i = i+1;
    end
end

%%
%{
for i = 1:length(results)
    pro = Problem{problem}('M',Objectives(obj),'D',(Objectives(obj)+10-1));
    cur_result = results{i,2};
    
    PopObj = cur_result.best.objs;
    fprintf("Sol: %d : %d %d \n", i, size(PopObj));
    isempty(PopObj)
    
    %gd  = GD(cur_result,pro.optimum);
    %igd  = IGD(cur_result,pro.optimum);
    %hv  = HV(cur_result,pro.optimum);
    %fprintf("gd %d, idg %d, hv %d \n",gd,igd,hv);
end
%}
%% remain on main directory
cd(my_dir);