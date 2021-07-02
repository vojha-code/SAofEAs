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

save("Testbench.mat", "problems_f");

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