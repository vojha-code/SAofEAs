
% f1 - Sphere Model
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -100, 'upper_bound', 100);
problem.obj_func = @(sol) f1(sol.x);
problems_f(1) = problem;
% f2 - Schwefel's Problem 2.22
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -10, 'upper_bound', 10);
problem.obj_func = @(sol) f2(sol.x);
problems_f(2) = problem;
% f3 - Schwefel's Problem 1.2
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -100, 'upper_bound', 100);
problem.obj_func = @(sol) f3(sol.x);
problems_f(3) = problem;
% f4 - Schwefel's Problem 2.21
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -100, 'upper_bound', 100);
problem.obj_func = @(sol) f4(sol.x);
problems_f(4) = problem;
% f5 - Generalized Rosenbrock Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -30, 'upper_bound', 30);
problem.obj_func = @(sol) f5(sol.x);
problems_f(5) = problem;
% f6 - Step Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -100, 'upper_bound', 100);
problem.obj_func = @(sol) f6(sol.x);
problems_f(6) = problem;
% f7 - Quartic Function w/ Noise
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -1.28, 'upper_bound', 1.28);
problem.obj_func = @(sol) f7(sol.x);
problems_f(7) = problem;
% f8 - Generalized Schwefel's Problem 2.26
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -500, 'upper_bound', 500);
problem.obj_func = @(sol) f8(sol.x);
problems_f(8) = problem;
% f9 - Generalized Rastrigrin Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -5.12, 'upper_bound', 5.12);
problem.obj_func = @(sol) f9(sol.x);
problems_f(9) = problem;
% f10 - Ackley's Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -32, 'upper_bound', 32);
problem.obj_func = @(sol) f10(sol.x);
problems_f(10) = problem;
% f11 - Generalized Griewank Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -600, 'upper_bound', 600);
problem.obj_func = @(sol) f11(sol.x);
problems_f(11) = problem;
% f12 - Generalized Penalized Function 1
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -50, 'upper_bound', 50);
problem.obj_func = @(sol) f12(sol.x);
problems_f(12) = problem;
% f13 - Generalized Penalized Function 2
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 30, 'lower_bound', -50, 'upper_bound', 50);
problem.obj_func = @(sol) f13(sol.x);
problems_f(13) = problem;
% f14 - Shekel's Foxholes function
holes = foxholes([-32,-16,0,16,32],5);
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 2, 'lower_bound', -65536, 'upper_bound', 65536);
problem.obj_func = @(sol) f14(sol.x, holes);
problems_f(14) = problem;
% f15 - Kowalik's function
kowalik_table.i = 1:11;
kowalik_table.a = [0.1957,0.1947,0.1735,0.1600,0.0844,0.0627,0.0456,0.0342,0.0323,0.0235,0.0246];
kowalik_table.b = 1./[0.25,0.5,1,2,4,6,8,10,12,14,16]; 
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 4, 'lower_bound', -5, 'upper_bound', 5);
problem.obj_func = @(sol) f15(sol.x, kowalik_table);
problems_f(15) = problem;
% f16 - Six-Hump Camel-Back Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 2, 'lower_bound', -5, 'upper_bound', 5);
problem.obj_func = @(sol) f16(sol.x);
problems_f(16) = problem;
% f17 - Branin Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 2, 'lower_bound', [-5,0], 'upper_bound', [10,15]);
problem.obj_func = @(sol) f17(sol.x);
problems_f(17) = problem;
% f18 - Goldstein-Price Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 2, 'lower_bound', -2, 'upper_bound', 2);
problem.obj_func = @(sol) f18(sol.x);
problems_f(18) = problem;
% f19 - Hartman Function (n=3)
hartman_table.a = [ ...
    3,10,30;
    0.1,10,35;
    3,10,30;
    0.1,10,35];
hartman_table.p = [ ...
    0.3689,0.1170,0.2673;
    0.4699,0.4387,0.7470;
    0.1091,0.8732,0.5547;
    0.038150,0.5743,0.8828];
hartman_table.c=[1,1.2,3,3.2];
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 3, 'lower_bound', 0, 'upper_bound', 1);
problem.obj_func = @(sol) hartman(sol.x, hartman_table);
problems_f(19) = problem;
% f20 - Hartman Function (n=6)
hartman_table.a = [ ...
    10,3,17,3.5,1.7,8;
    0.05,10,17,0.1,8,14;
    3,3.5,1.7,10,17,8;
    17,8,0.05,10,0.1,14];
hartman_table.p = [ ...
    0.1312,0.1696,0.5569,0.0124,0.8283,0.5886;
    0.2329,0.4135,0.8307,0.3736,0.1004,0.9991;
    0.2348,0.1415,0.3522,0.2883,0.3047,0.6650;
    0.4047,0.8828,0.8732,0.5743,0.1091,0.0381];
hartman_table.c=[1,1.2,3,3.2];
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 6, 'lower_bound', 0, 'upper_bound', 1);
problem.obj_func = @(sol) hartman(sol.x, hartman_table);
problems_f(20) = problem;
% f21 - Shekel Function(m=5)
shekel_table.a=[...
    4,4,4,4;
    1,1,1,1;
    8,8,8,8;
    6,6,6,6;
    3,7,3,7;
    2,9,2,9;
    5,5,3,3;
    8,1,8,1;
    6,2,6,2;
    7,3.6,7,3.6];
shekel_table.c=[ 0.1,0.2,0.2,0.4,0.4,0.6,0.3,0.7,0.5,0.5];
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 4, 'lower_bound', 0, 'upper_bound', 10);
problem.obj_func = @(sol) shekel(sol.x, shekel_table, 5);
problems_f(21) = problem;
% f22 - Shekel Function(m=7)
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 4, 'lower_bound', 0, 'upper_bound', 10);
problem.obj_func = @(sol) shekel(sol.x, shekel_table, 7);
problems_f(22) = problem;
% f23 - Shekel Function(m=10)
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 4, 'lower_bound', 0, 'upper_bound', 10);
problem.obj_func = @(sol) shekel(sol.x, shekel_table, 10);
problems_f(23) = problem;



% SR 1 Shifted Shpere:  sphere(M*(x-o))
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f24(sol.x);
problems_f(24) = problem;

% SR 2 Shifted Ellipsoid function, we had better not shift and rotate it 
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f25(sol.x);
problems_f(25) = problem;

% SR 3 Shifted Ackley, we shift it. Ackley(x-o)
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -32, 'upper_bound', 32);
problem.obj_func = @(sol) f26(sol.x);
problems_f(26) = problem;

% SR 4 Shifted Griewank, we shift it. Griewank(x-o)
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -50, 'upper_bound', 50);
problem.obj_func = @(sol) f27(sol.x);
problems_f(27) = problem;

% SR 5 Shifted and Rotated Rosenbrock, we shift and rotate it. Rosenbrock(M*(x-o))
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f28(sol.x);
problems_f(28) = problem;

% SR 6 Shifted and Rotated Rastrigin, we shift and rotate it. Rastrigin(M*(x-o))
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f29(sol.x);
problems_f(29) = problem;

% SR 7 Shifted and Rotated Weierstrass Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -100, 'upper_bound', 100);
problem.obj_func = @(sol) f30(sol.x);
problems_f(30) = problem;

% SR 8 Shifted and Rotated Schwefel’s Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f31(sol.x);
problems_f(31) = problem;

% SR 9 Shifted and Rotated Katsuura Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f32(sol.x);
problems_f(32) = problem;

% SR 10 Shifted and Rotated HappyCat Function
problem = ypea_problem();
problem.type = 'min';
problem.vars = ypea_var('x', 'real', 'size', 10, 'lower_bound', -20, 'upper_bound', 20);
problem.obj_func = @(sol) f33(sol.x);
problems_f(33) = problem;

save("Testbench.mat", "problems_f");
%% NEW funcitons from CEC
function z = f24(x)  
    %[row, col] = size(x);    
    %dim =  max(row,col);
    %x = reshape(x, dim, 1);
    sd =  [50.3557898229086,64.9267099320991,-59.6821093930390,66.1401369822431,21.1774793960655,-64.3935352000945,-35.4402849812723,7.50104307279742,73.2010936694876,74.3821656318843];
    %sd = load('sd_1_D10.mat').sd;
    %x_rs = x - reshape(sd(1:dim), size(x)) / 8;
    x_rs = x - sd / 8;
    %x = x_rs;
    z = sum(x_rs.*x_rs);
end

function z = f25(x)
    %[row, col] = size(x);    
    %x = reshape(x, max(row,col), 1);
    %dim = max(row,col);
    
    %sd = load('sd_2_D10.mat').sd;
    sd = [-54.0508306890812,47.0855265094251,-30.2055932728312,4.56530168099404,-53.4962032800351,16.3171106242618,-37.9245944735769,24.6526557562852,30.2743205024012,39.7042548517935];
    %x_rs = x - reshape(sd(1:dim), size(x)) / 8;
    
    x_rs = x - sd / 8;
    %x = x_rs;
    %z = sum(x .* x .* [1:numel(x)]');
    
    z = sum(x_rs .* x_rs * [1:numel(x_rs)]');
    %size(z)
end
function z = f26( x )
    %[row, col] = size(x);    
    %x = reshape(x, max(row,col), 1);
    %dim = max(row,col);
    
    %sd = load('sd_3_D10.mat').sd;
    sd = [-12.3382897439865,-64.9233057779625,15.7637870010786,-4.65211898266657,31.3518901282573,31.9820559885267,22.1649213234941,-74.6233862293713,-68.9910241411118,-28.8640423711206];
    %x_rs = x - reshape(sd(1:dim), size(x)) / 8;
    x_rs = x - sd / 8;
    
    %x = x_rs;
    
    %z = -20 * exp(-0.2 * sqrt( sum(x .* x)/ numel(x))) ...
    %    - exp(sum(cos(2*pi*x)) / numel(x)) ...
    %    + 20 + 2.7182818284590452353602874713526625;
    
    z = -20 * exp(-0.2 * sqrt( sum(x_rs .* x_rs)/ numel(x_rs))) ...
        - exp(sum(cos(2*pi*x_rs)) / numel(x_rs)) ...
        + 20 + 2.7182818284590452353602874713526625;
end

function z = f27( x )
    %[row, col] = size(x);    
    %x = reshape(x, max(row,col), 1);
    %dim = max(row,col);
    
    %sd = load('sd_4_D10.mat').sd;
    sd =  [13.1598663243563,6.51829393990555,59.1905651772812,-37.6353557638992,-29.1081479230305,-60.9256734313294,70.3727152551873,23.2882999956038,-3.28588400817795,22.2907137664173];
    %x_rs = x - reshape(sd(1:dim), size(x)) / 8;
    x_rs = x - sd / 8;
    %x = x_rs;
    
    %z = sum(x.^2) ./ 4000 - prod(cos(x ./ sqrt([1:numel(x)]'))) + 1;
    z = sum(x_rs.^2) ./ 4000 - prod(cos(x_rs ./ sqrt([1:numel(x_rs)]))) + 1;
end

function z = f28(x)
    %[row, col] = size(x);    
    %x = reshape(x, max(row,col), 1);
    %dim = max(row,col);
    
    %sd = load('sd_5_D10.mat').sd;
    sd = [56.1140278862412,9.68952437678156,68.7374186810661,31.4667520888364,13.2465544281344,50.4635538363874,60.6422247355484,78.2258585727343,-79.9164199428888,58.4701745620839];
    %ma = load('m_5_D10.mat').ma;
    ma = [0.418619762832516,-0.721758017164471,0,0,0,0,-0.551201287031207,0,0,0;0.471931370056398,-0.345668776995509,0,0,0,0,0.811044930054507,0,0,0;0,0,-0.601669102579926,0,1.83083762297913,0.291203961012283,0,0,0,0;0,0,0,-0.950465563497317,0,0,0,-0.229833311446323,0.0486840357867678,0.203523281751394;0,0,-1.20205824023904,0,-0.377612080345293,-0.451207268417661,0,0,0,0;0,0,-0.998044682836883,0,-0.708838101293890,0.947753096093688,0,0,0,0;-0.775911255313866,-0.599648614831420,0,0,0,0,0.195916468445152,0,0,0;0,0,0,-0.0187997580682718,0,0,0,0.248516598721124,-0.880407501297162,0.403446032222832;0,0,0,-0.192563162863654,0,0,0,-0.0748270757450348,-0.421611508735615,-0.882929256931194;0,0,0,0.243271884956417,0,0,0,-0.937985638400263,-0.211556687487961,0.127457836650892];
    %x_rs = x - reshape(sd(1:dim), size(x)) / 8;
    x_rs = x - sd / 8;
    x_rs = 2.08 / 20 * ma * x_rs' + 1;
    %x = x_rs;
    
    %z = sum(100*(x(1:end-1) .^2 - x(2:end)) .^ 2 + (x(1:end-1)-1) .^2);
    z = sum(100*(x_rs(1:end-1) .^2 - x_rs(2:end)) .^ 2 + (x_rs(1:end-1)-1) .^2);
end


function z = f29(x)
    %[row, col] = size(x);    
    %x = reshape(x, max(row,col), 1);
    %dim = max(row,col);
    
    %sd = load('sd_6_D10.mat').sd;
    sd = [9.44520719819121,56.6559918837509,-24.3393289076382,-8.63573631118355,-71.3216824894193,-51.6627945937425,26.0492899137559,-27.0673607674713,63.7577820534879,-61.0951682485263];
    %ma = load('m_6_D10.mat').ma;
    ma = [0.871816330973528,0,0,0,0,-1.32147073673012,0,0,0,0.409793616048068;0,0.484376173326003,-0.0378558565124036,0,-0.106672391979856,0,0,-0.867506574978579,0,0;0,0.384178794355382,0.402794164116801,0,-0.777549217839257,0,0,0.292541841722945,0,0;0,0,0,-0.745876201419315,0,0,-0.410165810707159,0,0.524816825076362,0;0,-0.720697078077999,0.517443509420397,0,-0.236968588933656,0,0,-0.395845707400409,0,0;0.481811115159795,0,0,0,0,0.108612877934832,0,0,0,1.73737728160317;0,0,0,-0.277382155892717,0,0,-0.525073205124831,0,-0.804585153263622,0;0,-0.313659207311913,-0.754039793406066,0,-0.572611650672293,0,0,-0.0718177494979577,0,0;0,0,0,0.605580574117600,0,0,-0.745695740223100,0,0.277866930858138,0;-0.875788902499223,0,0,0,0,-0.384951297393576,0,0,0,0.339617919552162];
    %x_rs = x - reshape(sd(1:dim), size(x)) / 8;
    x_rs = x - sd / 8;
    x_rs = 5.12 / 20 * ma * x_rs';
    %x = x_rs;
    %z = sum(x .^ 2 - 10 * cos(2*pi*x) + 10);
    z = sum(x_rs .^ 2 - 10 * cos(2*pi*x_rs) + 10);
end

function z =  f30(x)
    nx =  length(x);
    %sd = load('sd_7_D10.mat').sd;
    sd = [-18.6856040000000,-59.0481370000000,-29.7200450000000,-61.8632000000000,59.1839890000000,37.6977510000000,25.3661660000000,-55.6511790000000,-73.1508630000000,-51.9150100000000];
    %ma = load('m_7_D10.mat').ma;
    ma = [0.364824000000000,0.327116000000000,0.331429000000000,0.383080000000000,0.305119000000000,0.296899000000000,0.159168000000000,0.273838000000000,0.300114000000000,0.362858000000000;0.205968000000000,-0.316605000000000,0.778059000000000,-0.104924000000000,0.0500950000000000,-0.125177000000000,-0.117077000000000,-0.422545000000000,0.0786000000000000,-0.156033000000000;-0.172747000000000,-0.0198120000000000,0.0251820000000000,-0.309099000000000,-0.0430160000000000,0.455149000000000,0.737796000000000,-0.172768000000000,0.213751000000000,-0.211417000000000;-0.289777000000000,0.359191000000000,0.176203000000000,-0.0161800000000000,-0.512633000000000,0.298378000000000,-0.446456000000000,0.0315970000000000,0.420004000000000,-0.164787000000000;0.271750000000000,-0.515025000000000,-0.398445000000000,0.414867000000000,-0.0278850000000000,0.313776000000000,-0.159061000000000,-0.252887000000000,0.347843000000000,-0.143348000000000;0.0772610000000000,0.415795000000000,-0.229569000000000,0.00482700000000000,0.154217000000000,-0.540508000000000,0.116497000000000,-0.433952000000000,0.495817000000000,-0.0690430000000000;-0.162683000000000,0.324026000000000,0.0462340000000000,0.400401000000000,0.375644000000000,0.187492000000000,-0.0657990000000000,-0.173513000000000,-0.361791000000000,-0.603732000000000;0.503246000000000,0.0263120000000000,0.0293620000000000,0.0346030000000000,-0.346122000000000,-0.245942000000000,0.202733000000000,0.458144000000000,0.0278930000000000,-0.558506000000000;-0.526382000000000,-0.341965000000000,0.124739000000000,0.134200000000000,0.340880000000000,-0.250533000000000,0.0293800000000000,0.444349000000000,0.402878000000000,-0.181182000000000;-0.272124000000000,-0.0508390000000000,0.151915000000000,0.628308000000000,-0.487935000000000,-0.225786000000000,0.363492000000000,-0.160953000000000,-0.144398000000000,0.193836000000000];
    sh_rate =  0.5/100.0;
    xshift = (x - sd)*sh_rate;
    srx = (ma*xshift')*sh_rate;
    %x = srx;
    %nx = dim;
    
    a = 0.5;
    b = 3.0;
    k_max = 20;
    f0 = 0.0;
    
	for i =1:nx
        zi = srx(i);
		sum = 0.0;
		sum2 = 0.0;
		for j = 0:k_max
			sum = sum + (a^j)*cos(2.0*pi*(b^j)*(zi+0.5));
			sum2 = sum2 + (a^j)*cos(2.0*pi*(b^j)*0.5);
        end
		f0 = f0 + sum;
    end
	z = f0 - nx*sum2;
end

function z = f31(x) 
    d =  length(x);
    %sd = load('sd_8_D10.mat').sd;
    sd = [-23.2522970000000,-40.9347650000000,-54.3118540000000,32.9758060000000,37.1078810000000,-79.9037050000000,36.0005540000000,-37.4181970000000,-51.3373520000000,76.1581060000000];
    %ma = load('m_8_D10.mat').ma;
    ma =  [0.447843000000000,0.443329000000000,0.391889000000000,0.118429000000000,0.0739630000000000,0.191381000000000,0.0359140000000000,0.298078000000000,0.359390000000000,0.417007000000000;-0.0531420000000000,0.316462000000000,-0.120080000000000,0.636874000000000,0.00983500000000000,-0.159254000000000,0.357002000000000,0.353204000000000,-0.169804000000000,-0.412921000000000;0.101147000000000,-0.581170000000000,0.347855000000000,0.343909000000000,0.391490000000000,0.291554000000000,0.362542000000000,-0.206687000000000,-0.0143080000000000,0.0102610000000000;-0.173198000000000,0.165060000000000,-0.691848000000000,0.0502190000000000,0.253953000000000,0.313921000000000,0.306709000000000,-0.170532000000000,0.268721000000000,0.321216000000000;-0.620070000000000,0.0625260000000000,0.209266000000000,-0.0208890000000000,0.0267080000000000,0.518586000000000,-0.148556000000000,0.389464000000000,-0.309785000000000,0.167370000000000;0.0390110000000000,0.490708000000000,0.228871000000000,-0.308266000000000,0.0779100000000000,0.354628000000000,0.225997000000000,-0.457535000000000,-0.129282000000000,-0.448683000000000;-0.383683000000000,0.113090000000000,0.224817000000000,-0.274749000000000,0.605491000000000,-0.509565000000000,0.161079000000000,0.110688000000000,0.222472000000000,0.000318000000000000;-0.370767000000000,0.217848000000000,0.255052000000000,0.432049000000000,-0.260387000000000,-0.233847000000000,-0.0625970000000000,-0.563228000000000,0.0237760000000000,0.345198000000000;0.0178680000000000,0.0797300000000000,-0.0852050000000000,0.319150000000000,0.429947000000000,0.142451000000000,-0.729541000000000,-0.121605000000000,0.215785000000000,-0.292368000000000;0.292191000000000,0.164319000000000,-0.115405000000000,0.00979700000000000,0.388094000000000,-0.158545000000000,-0.102542000000000,-0.0964880000000000,-0.748523000000000,0.344014000000000];
    sh_rate =  1000.0/100.0;
    xshift = (x - sd)*sh_rate;
    srx = (ma*xshift')*sh_rate;
    %x = srx;
    
    f0 = 0.0;
    for i = 1:d
        zi = srx(i) + 418.9829;
        if(zi > 500)
            f0 = f0 - (500 - mod(zi,500))*sin(sqrt(abs(500 - mod(zi,500))));
            temp = (zi - 500)/100;
            f0 = f0 + temp*temp/d;
        elseif (zi <-500)
            f0 = f0 - (-500 + mod(abs(zi),500))*sin(sqrt(500 - mod(abs(zi),500)));
            temp = (zi + 500)/100;
            f0 = f0 + temp*temp/d;
        else
            f0 = f0 - zi*sin(sqrt(abs(zi)));
        end
    end
    z = f0 + 418.9829*d;
end

function z = f32(x) 
    nx =  length(x);
    %sd = load('sd_9_D10.mat').sd;
    sd = [75.4904150000000,49.1138380000000,41.8540690000000,-59.0793190000000,29.6608120000000,52.7651920000000,-61.7357430000000,37.5556320000000,75.2577790000000,-53.5098150000000];
    %ma = load('m_9_D10.mat').ma;
    ma = [0.341000000000000,0.731029000000000,0.317794000000000,0.682806000000000,-0.357768000000000,-0.283196000000000,0.450559000000000,0.179870000000000,0.794492000000000,-0.174687000000000;-0.128909000000000,0.0364280000000000,0.634555000000000,-0.194699000000000,-0.106677000000000,-0.0345890000000000,0.799824000000000,-0.122041000000000,-0.798925000000000,-0.0735260000000000;-0.346802000000000,0.329111000000000,0.0740290000000000,-0.193849000000000,-0.434083000000000,0.0652870000000000,0.371303000000000,1.05337400000000,0.549114000000000,0.532533000000000;-0.0737630000000000,0.490117000000000,-0.362274000000000,0.493554000000000,-0.626449000000000,-0.484650000000000,-0.241967000000000,0.0599630000000000,-0.719277000000000,0.506735000000000;0.0293780000000000,0.00525300000000000,-0.0232260000000000,-0.490425000000000,-0.161565000000000,-0.0436460000000000,0.107907000000000,-0.712226000000000,0.205402000000000,0.845720000000000;-0.353685000000000,0.168015000000000,-0.764824000000000,-0.726978000000000,-0.416118000000000,-0.577103000000000,0.475662000000000,-0.180706000000000,-0.0694590000000000,-0.628644000000000;0.445996000000000,1.03630900000000,-0.0347610000000000,-0.383024000000000,0.270063000000000,0.795308000000000,-0.235941000000000,0.219084000000000,-0.288884000000000,0.00924700000000000;0.174309000000000,0.0267100000000000,-0.404853000000000,0.598886000000000,0.774549000000000,0.643057000000000,0.680597000000000,0.168621000000000,0.0488740000000000,0.210643000000000;1.05486000000000,-0.242008000000000,-0.0179300000000000,-0.213661000000000,-0.404517000000000,0.104756000000000,0.599850000000000,0.411779000000000,-0.0234010000000000,0.386477000000000;0.252959000000000,0.221437000000000,0.167872000000000,-0.165260000000000,1.23758400000000,-0.561960000000000,0.212229000000000,0.184356000000000,-0.590932000000000,-0.0490600000000000];
    sh_rate =  0.5/100.0;
    xshift = (x - sd)*sh_rate;
    srx = (ma*xshift')*sh_rate;
    %x = srx;
    
    %nx = dim;
    f0 = 1.0;
    tmp3 = (1.0*nx)^1.2;
    for i =1:nx
        zi = srx(i);
		temp = 0.0;
		for j =1:32
			tmp1 = 2.0^j;
			tmp2 = tmp1*zi;
			temp = temp + abs(tmp2 - floor(tmp2 + 0.5))/tmp1;
        end
		f0 = f0 * (1.0 + i*temp)^(10.0 / tmp3);
    end
	tmp1 = 10.0/nx/nx;
    z = f0*tmp1 - tmp1;
end


function z = f33(x) 
    nx =  length(x);
    %sd = load('sd_10_D10.mat').sd;
    sd = [-20.4240720000000,-18.0470520000000,70.0681490000000,-59.6251430000000,-7.75597000000000,24.1473970000000,-39.3589870000000,53.1425080000000,-41.5851780000000,6.16907000000000];
    %ma = load('m_10_D10.mat').ma;
    ma = [0.0210340000000000,0.261584000000000,0.485651000000000,0.776701000000000,-0.708169000000000,0.577365000000000,-0.797967000000000,0.398355000000000,-0.156178000000000,-0.344946000000000;0.261020000000000,0.258659000000000,0.210691000000000,-0.813783000000000,0.267241000000000,-0.0670860000000000,-0.144027000000000,0.651012000000000,0.0712570000000000,1.07427300000000;-0.311144000000000,0.568546000000000,0.174843000000000,0.414431000000000,-0.0269990000000000,-0.521969000000000,0.162726000000000,0.691829000000000,0.629823000000000,-0.463148000000000;0.775411000000000,-0.595822000000000,0.659629000000000,0.0880390000000000,0.550237000000000,0.179937000000000,0.202120000000000,-0.184791000000000,0.225303000000000,-0.545002000000000;-0.394878000000000,-0.583562000000000,-0.611139000000000,0.714331000000000,0.272331000000000,0.688121000000000,0.134190000000000,0.0607350000000000,0.617454000000000,0.0400500000000000;0.0515080000000000,-0.733838000000000,0.159088000000000,-0.228221000000000,-0.551424000000000,-0.303141000000000,-0.175428000000000,0.545259000000000,0.306215000000000,-0.265324000000000;-0.675659000000000,0.465085000000000,-0.0470700000000000,-0.747703000000000,0.0146750000000000,0.707116000000000,-0.0577760000000000,-0.0418460000000000,0.523591000000000,-0.424228000000000;-0.261814000000000,-0.135723000000000,-0.0304670000000000,-0.0465300000000000,0.811453000000000,-0.572000000000000,-0.737485000000000,0.252620000000000,0.0384630000000000,-0.401143000000000;-0.378889000000000,-0.106336000000000,0.127229000000000,-0.294674000000000,0.390475000000000,0.0264190000000000,0.752270000000000,0.847978000000000,-0.894553000000000,-0.407500000000000;-0.752360000000000,-0.339830000000000,0.657006000000000,-0.00586000000000000,0.143468000000000,-0.277829000000000,-0.00567400000000000,-0.226788000000000,0.522272000000000,0.506995000000000];
    sh_rate =  0.5/100.0;
    xshift = (x - sd)*sh_rate;
    srx = (ma*xshift')*sh_rate;
    %x = srx;
    
    %nx = dim;    
    alpha = 1.0/8.0;
	r2 = 0.0;
	sum_z=0.0;
    for i = 1:nx
        zi = srx(i);
		zi = zi - 1.0;% shift to orgin
        r2 = r2 + zi*zi;
		sum_z = sum_z + zi;
    end
    z = abs(r2 - nx)^(2*alpha) + (0.5*r2 + sum_z)/nx + 0.5;
end

% functions
function z=f1(x)
    z = sum(x.^2);
end
function z=f2(x)
    z = sum(abs(x))+prod(abs(x));
end
function z=f3(x)
    z = sum(cumsum(x).^2);
end
function z=f4(x)
    z = max(abs(x));
end
function z=f5(x)
    z = sum(100*(x(2:end)-x(1:end-1).^2).^2+(x(1:end-1)-1).^2);
end
function z=f6(x)
    z = sum(floor(x+0.5).^2);
end
function z=f7(x)
    i = 1:numel(x);
    z = sum(i.*x+rand(size(x)));
end
function z=f8(x)
    z = -sum(x.*sin(sqrt(abs(x))));
end
function z=f9(x)
    z = sum(x.^2-10*cos(2*pi*x)+10);
end
function z=f10(x)
    z = 20*(1-exp(-0.2*sqrt(mean(x.^2))))+exp(1)-exp(mean(cos(2*pi*x)));
end
function z=f11(x)
    i = 1:numel(x);
    z = 1/4000 * sum(x.^2)-prod(cos(x/sqrt(i)))+1;
end
function z=f12(x)
   n=numel(x);
   y=1+1/4*(x+1);
   z=pi/n*( 10*sin(pi*y(1))^2 + ...
       sum((y(1:end-1)-1).^2.*(1+10*sin(pi*y(2:end)).^2))+(y(end)-1)^2) + ...
       sum(u(x,10,100,4));
end
function z=u(x,a,k,m)
    z=(abs(x)>a).*k.*(sign(x)-a).^m;
end
function z=f13(x)
    z=0.1*(sin(3*pi*x(1))^2 + ...
        sum((x(1:end-1)-1).^2.*(1+sin(3*pi*x(2:end))))) + ...
        (x(end)-1)*(1+sin(2*pi*x(end))^2)+sum(u(x,5,100,4));
end
function x=foxholes(array,repeat)
    x=zeros(2,repeat*numel(array));
    x(1,:)=repmat(array,1,repeat);
    x(2,:)=repelem(array,repeat);
end
function z=f14(x,a)
    j=1:size(a,2);
    z=(1/500 + sum((j+sum((x'-a).^6)).^(-1))).^(-1);
end
function z=f15(x,kowalik_table)
    z=sum((kowalik_table.a-...
        (x(1).*(kowalik_table.b.^2+kowalik_table.b.*x(2)))./...
        (kowalik_table.b.^2+kowalik_table.b.*x(3)+x(4))).^2);
end
function z=f16(x)
    z=4*x(1)^2-2.1*x(1)^4+1/3*x(1)^6+x(1)*x(2)-4*x(2)^2+4*x(2)^4;
end
function z=f17(x)
    z=(x(2)-5.1/(4*pi^2)*x(1)^2+5/pi*x(1)-6)^2 ...
    +10*(1-1/(8*pi))*cos(x(1))+10;
end
function z=f18(x)
    z=(1 + (x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2))*...
        (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2));
end
function z=hartman(x,hartman_table)
    z=-sum(hartman_table.c'.*exp(-sum(hartman_table.a.*(x-hartman_table.p).^2,2)));
end
function z=shekel(x,shekel_table,m)
    z=-sum((sum((x-shekel_table.a(1:m,:)).^2,2)'+shekel_table.c(1:m)).^(-1));
end