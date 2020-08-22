function [U1, U2, V1, V2, evalFull] = solve2Graphs(A1, A2, C0, param, d, e)

% Function that solves the system of differential equations determined by
% a reaction diffusion system following a Schnakenberg model on two weakly 
% connected graphs. 

% INPUT:
%    [A1] = adjacency matrix of one graph.
%    [A2] = adjacency matrix of the other graph.
%    [C0] = matrix of connections between the two graphs
%    [param] = 1 x 6 vector with the parameters of the Schnakenberg model 
%        in the following order: alpha1, beta1, zeta1, alpha2, beta2, zeta2.
%    [d] = diffusion coefficients du1 and dv1, du2, dv2.
%    [e] = value of epsilon.


% OUTPUT:
%  [U1] = N1 x 1 vector with the values of species u in each node of the
%  first graph.
%  [U2] = N2 x 1 vector with the values of species u in each node of the
%  second graph.
%  [V1] = N1 x 1 vector with the values of species v in each node of the
%  first graph.
%  [V2] = N2 x 1 vector with the values of species v in each node of the
%  second graph.

rng(100);       % set seed for reproducibility. 
N1 = length(A1); 
N2 = length(A2);
 
L1 = A1-diag(sum(A1,1));
L2 = A2-diag(sum(A2,1));

D1 = -diag(sum(C0,2));
D2 = -diag(sum(C0,1));

alpha1 = param(1);
beta1 = param(2);
zeta1 = param(3);

alpha2 = param(4);
beta2 = param(5);
zeta2 = param(6);

J1 = [-(beta1+zeta1)^2/alpha1^2, -2*beta1*alpha1/(beta1+zeta1);(beta1+zeta1)^2/alpha1^2, 2*beta1*alpha1/(beta1+zeta1)-alpha1];
J2 = [-(beta2+zeta2)^2/alpha2^2, -2*beta2*alpha2/(beta2+zeta2);(beta2+zeta2)^2/alpha2^2, 2*beta2*alpha2/(beta2+zeta2)-alpha2];

[evalFull, evecFull, evalAsym0, evecAsym0, evalAsym1, lambdaAsym] = ComputeGraphAsymptoticsSystem(A1, A2, C0, e, J1, J2, d);


IC1 = [beta1*alpha1^2./(beta1+zeta1).^2*ones(N1,1);(beta1+zeta1)./alpha1*ones(N1,1)].*(1+normrnd(0,1e-3,[2*N1,1]));  
f1 = @(U)beta1-U(N1+1:2*N1).^2.*U(1:N1); 
g1 = @(U)U(N1+1:2*N1).^2.*U(1:N1)-alpha1.*U(N1+1:2*N1)+zeta1;

IC2 = [beta2*alpha2^2./(beta2+zeta2).^2*ones(N2,1);(beta2+zeta2)./alpha2*ones(N2,1)].*(1+normrnd(0,1e-3,[2*N2,1]));  
f2 = @(U)beta2-U(2*N1+N2+1:end).^2.*U(2*N1+1:2*N1+N2); 
g2 = @(U)U(2*N1+N2+1:end).^2.*U(2*N1+1:2*N1+N2)-alpha2.*U(2*N1+N2+1:end)+zeta2;

tspan = linspace(0,1e6,10); 
myfun = @(t,U)[f1(U)+d(1)*L1*U(1:N1)+e*d(1)*(D1*U(1:N1)+C0*U(2*N1+1:2*N1+N2)); g1(U)+d(2)*L1*U(N1+1:2*N1)+e*d(2)*(D1*U(N1+1:2*N1)+C0*U(2*N1+N2+1:end));
    f2(U)+d(3)*L2*U(2*N1+1:2*N1+N2)+e*d(3)*(D2*U(2*N1+1:2*N1+N2)+C0'*U(1:N1)); g2(U)+d(4)*L2*U(2*N1+N2+1:end)+e*d(4)*(D2*U(2*N1+N2+1:end)+C0'*U(N1+1:2*N1))];

options = odeset('RelTol',1e-9,'AbsTol',1e-9); 

[t y] = ode15s(myfun, tspan, [IC1;IC2],options);    %solve the system of ODEs numerically. 
U1 = y(:,1:N1); 
V1 = y(:,N1+1:2*N1);
U2 = y(:,2*N1+1:2*N1+N2); 
V2 = y(:,2*N1+N2+1:end);
end