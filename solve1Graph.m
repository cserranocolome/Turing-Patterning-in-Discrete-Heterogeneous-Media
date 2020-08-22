function [U1,V1] = solve1Graph(A1, param, d)

% Function that solves the system of differential equations determined by
% a reaction diffusion system following a Schnakenberg model with adjacency
% matrix "A1" and kinetic parameters "param".

% INPUT:
%    [A1] = adjacency matrix of the graph.
%    [param] = 1 x 3 vector with the parameters of the Schnakenberg model in the following
%    order: alpha, beta, zeta.
%    [d] = diffusion coefficients d_1 and d_2.


% OUTPUT:
%  [U1] = N x 1 vector with the values of species u in each node of the
%  graph.
%  [V1] = N x 1 vector with the values of species v in each node of the
%  graph.

rng(100);   % set seed for reproducibility. 
N1 = length(A1); 
 
L1 = A1-diag(sum(A1,1));

alpha1 = param(1);
beta1 = param(2);
zeta1 = param(3);

IC1 = [beta1*alpha1^2./(beta1+zeta1).^2*ones(N1,1);(beta1+zeta1)./alpha1*ones(N1,1)].*(1+normrnd(0,1e-3,[2*N1,1]));  
f1 = @(U)beta1-U(N1+1:2*N1).^2.*U(1:N1); 
g1 = @(U)U(N1+1:2*N1).^2.*U(1:N1)-alpha1.*U(N1+1:2*N1)+zeta1;


tspan = linspace(0,1e6,10); 
myfun = @(t,U)[f1(U)+d(1)*L1*U(1:N1); g1(U)+d(2)*L1*U(N1+1:2*N1)];
   
options = odeset('RelTol',1e-9,'AbsTol',1e-9); 

[t y] = ode15s(myfun, tspan, IC1,options);  %solve the system of ODEs numerically. 
U1 = y(:,1:N1); 
V1 = y(:,N1+1:2*N1);

end