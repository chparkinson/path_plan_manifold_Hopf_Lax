clear; %close all;

% for reproducible results
SEED = 501;
rng(SEED);

%choose horizon time and discretization level
t = 4.5; J = round(10*t);

% choose how many paths you want (randomly generate initial points)
HowManyPaths = 10;
for i = 1:HowManyPaths
    x_target{i} = 2*rand(10,1)-1;
end

% choose final point
xf = zeros(10,1);

% set PDHG parameters
sig = 1; tau = 0.25/((1+2*pi^2)*sig); kappa = 1; max_iter = 40000; tol = 1e-3;

% run Algorithm 1 to generate each optimal path (record CPU time for each)
for i = 1:HowManyPaths
    fprintf("==================== Trial %02i =======================\n",i);
    TIMERRR(i) = tic;
    %%%%
    %%%% Here is where optimal paths are resolved
    %%%%
    [u(i),x{i},p{i},howManyIter(i)] = HJBSolve(x_target{i},xf,t,J,sig,tau,kappa,max_iter,tol);
    %%%%
    %%%%
    %%%%
    TIME(i) = toc(TIMERRR(i));
    if howManyIter(i) == max_iter
        fprintf("Failed to converge in %i iterations\n",max_iter);
    else
        fprintf("Pathfinder converged in %i iterations. CPU time: %.2f sec\n",howManyIter(i),TIME(i));
    end
end
fprintf("=====================================================\n",i);
%% report results
N = zeros(1,HowManyPaths);
for j = 1:HowManyPaths
    N(j) = norm(x_target{j},2);
end


% print table for LaTeX tabular environment
fprintf('\\abs{x} & u & Err. \\\\ \n') 
fprintf('\\hline\n');
for j = 1:10
   fprintf('%.4f  &  %.4f  & %.4e \\\\ \n', N(j),u(j),abs(N(j)-u(j))); 
end
fprintf('\\hline\n');

