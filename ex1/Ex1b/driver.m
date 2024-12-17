clear; %close all;

t = 4.5; J = round(10*t);

SEED = 249;
rng(SEED)

%%
a = 1;

x_target = [-1; -1];
xf = [1; 1];

sig = 1; tau = 0.25/((1+2*a^2*pi^2)*sig); theta = 1; max_iter = 50000; tol = 1e-3;
TIME = tic;
[u,x,p,howManyIter] = HJBSolve(x_target,xf,t,J,sig,tau,theta,max_iter,tol,a);
TIME = toc(TIME);
if howManyIter == max_iter
    fprintf("Failed to converge in %i iterations\n",max_iter);
else
    fprintf("Pathfinder converged in %i iterations. CPU time: %.2f sec\n",howManyIter,TIME);
end
%%
[X,Y] = ndgrid(-1.2:0.05:1.2,-1.2:0.05:1.2);
F = figure(21); clf; hold on; 
M = @(x,y) a*sin(pi*x).*cos(pi*y); 
surf(X,Y,M(X,Y),'edgecolor','none');
plot3(x(1,:),x(2,:),M(x(1,:),x(2,:))+0.05,'k','linewidth',2);
plot3(xf(1),xf(2),M(xf(1),xf(2))+0.05,'r.','markersize',20);
plot3(x_target(1),x_target(2),M(x_target(1),x_target(2))+0.05,'.','markersize',20,'color','g');
axis([-1.1 1.1 -1.1 1.1 -3.1 3.1]);
view([0 55]);
axis off;

%%% print picture if desired
% print('pic2','-dpng');
% clearvars F; 

%%% save results if desired
% save Ex1b1.mat;
%%

a = 3;

x_target = [-1; -1];
xf = [1; 1];

sig = 1; tau = 0.25/((1+2*a^2*pi^2)*sig); theta = 1; max_iter = 50000; tol = 1e-3;
TIME = tic;
[u,x,p,howManyIter] = HJBSolve(x_target,xf,t,J,sig,tau,theta,max_iter,tol,a);
TIME = toc(TIME);
if howManyIter == max_iter
    fprintf("Failed to converge in %i iterations\n",max_iter);
else
    fprintf("Pathfinder converged in %i iterations. CPU time: %.2f sec\n",howManyIter,TIME);
end
%%
[X,Y] = ndgrid(-1.2:0.05:1.2,-1.2:0.05:1.2);
F = figure(22); clf; hold on; 
M = @(x,y) a*sin(pi*x).*cos(pi*y); 
surf(X,Y,M(X,Y),'edgecolor','none');
plot3(x(1,:),x(2,:)+0.0165,M(x(1,:),x(2,:))+0.05,'k','linewidth',2);
plot3(xf(1),xf(2),M(xf(1),xf(2))+0.05,'r.','markersize',20);
plot3(x_target(1),x_target(2),M(x_target(1),x_target(2))+0.05,'.','markersize',20,'color','g');
axis([-1.1 1.1 -1.1 1.1 -3.1 3.1]);
view([0 55]);
axis off;

%%% print picture if desired
% print('pic3','-dpng');
% clearvars F; 

%%% save results if desired
% save Ex1b2.mat;