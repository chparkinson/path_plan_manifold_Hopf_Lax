function [u,x,p,howManyIter] = HJBSolve(x_target,xf,t,J,sig,tau,kappa,max_iter,tol)
% INPUT:
%   x_target - the x value where you will solve the HJB equation
%   xf - the endpoint that the path will travel toward
%   t - the time you will solve the HJB equation at
%   J - the discretization level
%   sig,tau,kappa - PDHG parameters
%   max_iter - maximum allowable iteration count
%   tol - convergence tolerance
% OUTPUT:
%   u - the value u(x_target,t) of the value function at the desired point
%   x - the optimal path in state space
%   p - the optimal path in co-state space
%   howManyIter - the number of iteration required for convergence
%                   (due to random initialization, this could change 
%                    somewhat significantly between runs)


% gradient descent parameters
gd_steps = 1;
gd_rate = 0.025;


% INITIALIZATION
dt = t/J;
dim = size(x_target,1);
x = zeros(dim,J+1);
for n = 1:J+1
    %initialize straight line path + randomness
    x(:,n) = (1-((n-1)/J))*xf + ((n-1)/J)*x_target + 0.2*randn(dim,1);
end
x(:,J+1) = x_target;
p = 0.1*randn(dim,J+1); p(:,1) = 0;
z = x;

% declaring variables 
GM = zeros(dim,J+1); % gradient of M along the path
HM = zeros(dim,dim,J+1); % hessian of M along the path
A = cell(1,J+1); % matrix A
V = zeros(1,J+1); % Velocity along the path
GV = zeros(dim,J+1); % gradient of velocity along the path
chi = zeros(1,J+1); % approx to indicator function
Gchi = zeros(dim,J+1); % gradient of approx to indicator function
B = 50; % parameter in approximation of indicator function
w=p;
for k = 1:max_iter
    pold = p; xold = x;

    for j = 2:J+1
        % here is where the manifold and velocity come in
        % VELOCITY
        V(j) = 1;
        GV(:,j) = [0;0];

        %% EXAMPLE 1(a): M = sin(pi*x)*cos(pi*y)
        GM(:,j) = pi*[cos(pi*x(1,j))*cos(pi*x(2,j));-sin(pi*x(1,j))*sin(pi*x(2,j))]; 
        HM(:,:,j) = -pi^2*[sin(pi*x(1,j))*cos(pi*x(2,j)),cos(pi*x(1,j))*sin(pi*x(2,j));
            cos(pi*x(1,j))*sin(pi*x(2,j)),sin(pi*x(1,j))*cos(pi*x(2,j))];

        %% Build matrices and approx to indicator
        A{j} = eye(dim) - GM(:,j)*GM(:,j)'/(1+GM(:,j)'*GM(:,j));
        chi(j) = 1-exp(-B*norm(x(:,j)-xf)^2);
        Gchi(:,j) = 2*B*exp(-B*norm(x(:,j)-xf)^2)*(x(:,j)-xf);
        % in the 2D examples, we hard code the (inverse) of the
        % Cholesky factorization
        Linv = (1/V(j))*[sqrt((1+norm(GM(:,j))^2)/(1+GM(2,j)^2)), 0;
            prod(GM(:,j))/sqrt(1+GM(2,j)^2), sqrt(1+GM(2,j)^2)]; 
        
        % Resolve updated p values
        b = w(:,j) + sig*(Linv*(z(:,j) - z(:,j-1)));
        w(:,j) = max(0,1-sig*dt*chi(j)*V(j)/norm(b,2))*b;
        p(:,j) = Linv'*w(:,j);
    end

    % Resolve updated x values
    if k <=2000
        % for the first 2000 iterations, we set x = nu
        x(:,1) = xf;
        x(:,2:J) = x(:,2:J)-tau*(p(:,2:J)-p(:,3:J+1));
        x(:,J+1) = x_target;
    else
        % after 2000 iterations, we perform gradient descent
        x(:,1) = xf;
        for j = 2:J
            nu = x(:,j)-tau*(p(:,j)-p(:,j+1));
            xs = x(:,j);
            for l = 1:gd_steps
                NORM = sqrt(p(:,j)'*A{j}*p(:,j));
                opgm2 = 1+GM(:,j)'*GM(:,j);
                IP = p(:,j)'*GM(:,j);
                GRAD = -dt*(Gchi(:,j)*(NORM-1)+(chi(j)/(NORM*opgm2^2+1e-8)*(IP^2*HM(:,:,j)*GM(:,j) - IP*opgm2*HM(:,:,j)*p(:,j)))) + (1/tau)*(xs-nu);
                xs = xs - gd_rate*GRAD;
            end
            x(:,j) = xs;
        end
        x(:,J+1) = x_target;
    end
    
    % update z values
    z = x + kappa*(x-xold);

    % compute change between iterations and break if small enough
    change = max(max(max(abs(x-xold))),max(max(abs(p-pold))));
    if change < tol
        break
    end

    % print progress 
    if mod(k,1000)==0
        fprintf('Iteration %i complete. Change = %.4e\n',k,change);
        % lower gradient descent rate and increase parameter B
        if k>=3000
            gd_rate = 0.5*gd_rate;
            if B<1000
                B = B+50;
            end
        end
    end
end

% compute u 
u = 0;
for j = 2:J+1
   V(j) = 1;
    GV(:,j) = [0;0];
    A{j} = eye(dim) - GM(:,j)*GM(:,j)'/(1+GM(:,j)'*GM(:,j));
    chi(j) = 1-exp(-B*norm(x(:,j)-xf)^2);
    Gchi(:,j) = 2*B*exp(-B*norm(x(:,j)-xf)^2)*(x(:,j)-xf);
    u = u + p(:,j)'*(x(:,j)-x(:,j-1))-dt*chi(j)*(V(j)*sqrt(p(:,j)'*A{j}*p(:,j))-1);
end

% record final iteration count
howManyIter = k;
end

