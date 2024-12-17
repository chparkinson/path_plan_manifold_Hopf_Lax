function [y] = plotHelp(x,EXTRA)
%A Helper function which takes a sparse path and 
% makes it finer using linear interpolation
% (this helps with plotting when the manifold is 
% concave and velocity is high; otherwise parts 
% of the path will be plotted "under" the manifold).

y = x(:,1);
for j = 2:size(x,2)
    x1 = x(:,j-1); x2 = x(:,j);
    for n = 1:EXTRA
        y(:,end+1) = x1 + (n/EXTRA)*(x2-x1);
    end
end
end

