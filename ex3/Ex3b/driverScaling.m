clear; %close all;
% This script runs scaling test for the high dimensional example

% fix rng to get results from paper (seed was chosen arbitrarily)
rng(9365);

% choose the dimensions you want to run the example in and 
% the number or trials per dimension
DIMS = 10:30;
TRIALS = 10;
TIME = zeros(TRIALS,length(DIMS));

% for each dimension, run the specified amount of trials and 
% record times
count = 1;
for dim = DIMS

    % make sure the time horizon is large enough based on dimension
    t = ceil(1.2*2*sqrt(dim)); N = round(10*t);

    x_target = -ones(dim,1); x_target(1) = -0.9;
    xf = ones(size(x_target));

    % run the code for the specified amount of trials and record times
    sig = 1; tau = 0.25/(10*sig); theta = 1; max_iter = 40000; tol = 1e-3;
    for m = 1:TRIALS
        fprintf("========= Dim = %i, Trial = %i =========\n",dim,m)
        TIMER = tic;
        [u,x,p,howManyIter] = HJBSolve(x_target,xf,t,N,sig,tau,theta,max_iter,tol);
        TIME(m,count) = toc(TIMER);
        if howManyIter == max_iter
            fprintf("Failed to converge in %i iterations\n",max_iter);
        else
            fprintf("Pathfinder converged in %i iterations. CPU time: %.2f sec\n",howManyIter,TIME(m,count));
        end
        fprintf("=======================================\n");
    end
    count = count+1;
end

% calc average and std of times
avgTIMES = mean(TIME);
stdevTIMES = std(TIME);

%% plot results
F = figure(1213);clf; hold on;
% plot(DIMS,avgTIMES,'k-','linewidth',2);
errorbar(avgTIMES,stdevTIMES,'linewidth',1);
axis([0 22 0.5 3.5])
xticks(1:4:21);
xticklabels(DIMS(1:4:end));
xlabel('Dimension');
ylabel('CPU time (sec)');
ax = gca; 
ax.FontSize = 20; 
ax.TickLabelInterpreter = 'latex';
ax.XLabel.Interpreter = 'latex';
ax.YLabel.Interpreter = 'latex';

%%% print picture if desired
% print('pic7','-dpng');

%%% save results if desired
% clearvars F;
% save Ex3b.mat