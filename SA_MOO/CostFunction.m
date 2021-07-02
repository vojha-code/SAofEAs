function f = CostFunction(x, nProblem)
    f = inf;
    if strcmp(nProblem,'MOP2')
        f = MOP2(x);  % Cost Function
    end
    %CostFunction = @(x) MOP2(x);  % Cost Function
end