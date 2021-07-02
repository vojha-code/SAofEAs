warning('off');

ProblemsList = { 
    %@ZDT1,
    %@ZDT2,
    %@ZDT3,
    %@ZDT4,
    %@ZDT5,
    %@ZDT6,
    @DTLZ1,
    @DTLZ2,
    @DTLZ3,
    @DTLZ4,
    %@DTLZ5,
    %@DTLZ6,
    %@DTLZ7,
    %@DTLZ8,
    %@DTLZ9,
    @C1_DTLZ1,
    %@C1_DTLZ3,
    %@C2_DTLZ2,
    %@C3_DTLZ4,
    @CDTLZ2,
    %@DC1_DTLZ1,
    %@DC1_DTLZ3,
    %@DC2_DTLZ1,
    %@DC2_DTLZ3,
    %@DC3_DTLZ1,
    %@DC3_DTLZ3,
    @IDTLZ1,
    @IDTLZ2,
    %@SDTLZ1,
    %@SDTLZ2,
    %@WFG1,
    %@WFG2,
    @WFG3,
    %@WFG4,
    %@WFG5,
    @WFG6,
    @WFG7,
    %@WFG8,
    %@WFG9
    };


ProblemsList = { 
    % time consuming problems
    @DTLZ2,
    @DTLZ4,
    @IDTLZ2,
    @WFG7,
    % easier problems
    @DTLZ1,
    @DTLZ3,
    @C1_DTLZ1,
    @CDTLZ2,
    @IDTLZ1,
    @WFG3,
    @WFG6
    };

%{
Dimensions = [3, 5, 10, 15, 8];  % Number of objectives
Problems = [];
for dim =1:length(Dimensions)
    for index =1:length(ProblemsList)
        Problems = [Problems; {ProblemsList{index},Dimensions(dim)}];
    end
end

%{
for i =1:length(Problems)
    Problem={Problems{i,1},Problems{i,2}};
    fprintf(" %s (%d) \n", func2str(Problem{1,1}),Problem{1,2});
end
%}

%}