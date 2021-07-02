
N = 100;
M = 2;
D = 10;
save = 0;
ITR = 1000;
maxFE = N*ITR;

platemo('algorithm',@MOEAD,'problem',@ZDT1,'N',N,'M',M,'D',D,'save',save,'maxFE',maxFE); 
%platemo('algorithm',@MOEAD,'problem',@ZDT1,'N',100,'M',2,'D',10,'save',0,'maxFE',(100*1000))
% pro = ZDT1('M',2,'D',10);
% igd  = IGD(result{end},pro.optimum);
% gd  = GD(result{end},pro.optimum);
% hv  = HV(result{end},pro.optimum);
%%
for i = 1:length(result)   
    pro = ZDT1('M',M,'D',D);
    cur_result = result{i,2};
    gd  = GD(cur_result,pro.optimum);
    igd  = IGD(cur_result,pro.optimum); 
    hv  = HV(cur_result,pro.optimum);
    fprintf("gd %d, idg %d, hv %d \n",gd,igd,hv);
end