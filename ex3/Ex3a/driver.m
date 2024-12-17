clear; %close all;
% This script runs the code for the high dimensional example

%fix rng for results from the paper (seed was chosen arbitrarily)
rng(2123);

% choose dimension and travel time / discretization level
dim = 25; 
t = ceil(1.2*2*sqrt(dim)); N = round(10*t);

% choose starting and ending point
x_target = -ones(dim,1); x_target(1) = -0.9; 
xf = ones(size(x_target));

% fix PDHG parameters
sig = 1; tau = 0.25/(10*sig); theta = 1; max_iter = 40000; tol = 1e-3;

% solve HJB equation and record time
TIME = tic;
[u,x,p,howManyIter] = HJBSolve(x_target,xf,t,N,sig,tau,theta,max_iter,tol);
TIME = toc(TIME);
if howManyIter == max_iter
    fprintf("Failed to converge in %i iterations\n",max_iter);
else
    fprintf("Pathfinder converged in %i iterations. CPU time: %.2f sec\n",howManyIter,TIME);
end
%% plot results
[X,Y] = ndgrid(-1.5:0.01:1.5,-1.5:0.01:1.5);
F = figure(21); clf; hold on;
plot(linspace(0,t,N+1),fliplr(x(1,:)),'color','k','LineWidth',2);
color = 0.5+0.5*rand(3,1); color(randi(3))=0;
plot(linspace(0,t,N+1),fliplr(x(2,:)),'color',color,'LineWidth',2);
for i = 3:length(x_target)
    color = 0.5+0.5*rand(3,1); color(randi(3))=0;
    plot(linspace(0,t,N+1),fliplr(x(i,:)),'color',color,'LineWidth',2,'HandleVisibility','off');
end
plot(linspace(0,t,N+1),fliplr(x(1,:)),'color','k','LineWidth',2,'HandleVisibility','off');
plot([u u],[-1.2,1.2],'--','color',[0.7 0.7 0.7],'linewidth',2,'HandleVisibility','off');
L = legend({'$x_1(t)$','$x_j(t)$'}); L.FontSize = 15; L.Interpreter ='latex'; L.Location ='northwest';
xlabel('$t$');
ax = gca; 
ax.FontSize = 20; 
ax.XLabel.Interpreter = 'latex';
ax.TickLabelInterpreter = 'latex';
axis([0 t -1.2 1.2]);

%%% print picture if desired
% print('pic6','-dpng');

%%% save results if desired
% clearvars F;
% save Ex3a.mat;
