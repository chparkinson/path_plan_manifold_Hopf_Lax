clear; %close all;

% for reproducible results
SEED = 501;
rng(SEED);

%choose horizon time and discretization level
t = 4.5; J = round(10*t);

% choose how many paths you want (randomly generate initial points)
HowManyPaths = 20;
for i = 1:HowManyPaths
    x_target{i} = 2*rand(2,1)-1;
end

% choose final point
xf = [1; 1];

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
%% plot results 
[X,Y] = ndgrid(-1.2:0.05:1.2,-1.2:0.05:1.2);
F = figure(21); clf; hold on;
M = @(x,y) sin(pi*x).*cos(pi*y);
surf(X,Y,M(X,Y),'edgecolor','none');
% for n = round(linspace(size(x{1},2),1,20))
for n = 1
    clf; hold on;
    surf(X,Y,M(X,Y),'edgecolor','none');
    for i = 1:HowManyPaths
        COLOR = 'k';
        % COLOR = 0.4*rand(3,1)+0.4; COLOR(randi(3))=0; % uncomment for more colorful paths :)
        plot3(x{i}(1,n:end),x{i}(2,n:end),M(x{i}(1,n:end),x{i}(2,n:end))+0.02,'color',COLOR,'linewidth',2);
        plot3(x{i}(1,end),x{i}(2,end),M(x{i}(1,end),x{i}(2,end))+0.05,'g.','markersize',20);
    end
    plot3(xf(1),xf(2),M(xf(1),xf(2))+0.03,'r.','MarkerSize',20);
    axis([-1.1 1.1 -1.1 1.1 -3.1 3.1]);
    view([0,55]);
    axis off;
    pause(0.1);
end

% Print picture if desired
% print('pic1','-dpng');

% save results if desired
% clearvars F;
% save Ex1a.mat;
