function [evalFull, evecFull, evalAsym0, evecAsym0, evalAsym1, lambdaAsym] = ComputeGraphAsymptotics(A1, A2, C0, e, fu1, fu2)

% Function that returns the asymptotic approximations of the eigenvalues of
% the full system composed by two coupled graphs. Scalar case.

% INPUT:
%    [A1] = N1 x N1 matrix. Adjacency matrix of the first graph.
%    [A2] = N2 x N2 matrix. Adjacency matrix of the second graph.
%  [C0] = N1 x N2 matrix. Coupling matrix
%   [e] = epsilon (strength of the coupling)
%   [fu1] = scalar f_u for the first graph
%   [fu2] = scalar f_u for the second graph

% OUTPUT:
%  [evalFull] = eigenvalues of the full system
%  [evecFull] = eigenvectors of the full system
%  [evalAsym0] = O(1) approximation of the eigenvalues
%  [evecAsym0] = O(1) approximation of the eigenvectors
%  [evalAsym1] = O(epsilon) term of the eigenvalues
%  [lambdaAsym] = O(epsilon) approximation of the eigenvalues:
%                   evalAsym0+epsilon*evalAsym1

N1 = size(A1,1);
N2 = size(A2,1);

D1 = -diag(sum(C0,2));
D2 = -diag(sum(C0,1));
D = zeros(N1+N2); D(1:N1,1:N1) = D1; D(N1+1:end,N1+1:end) = D2;

%Full Graph
C = e*C0;
A0 = [A1,0*C;0*C',A2];
A3 = [A1,C;C',A2];

%Full Laplacian
L0 = A0-diag(sum(A0,1));
L3 = A3-diag(sum(A3,1));

%Full system
M1 = fu1*eye(N1); M2 = fu2*eye(N2);
M = zeros(N1+N2); M(1:N1,1:N1) = M1; M(N1+1:end,N1+1:end) = M2;
M0 = M + L0;    %J+L
M = M+L3;       %J+L+extra connections

%Eigenvalues and eigenvectors of the full system
[VFull,evalFull] = eig(M);
[evalFull,indFull] = esort(diag(evalFull));
evecFull = VFull(:,indFull);

% Asymptotic evalues:
[V,evalAsym0] = eig(M0);
[evalAsym0,ind] = esort(diag(evalAsym0));
evecAsym0 = V(:,ind);
evalAsym1 = diag(evecAsym0'*D*evecAsym0);

lambdaAsym = evalAsym0 + e*evalAsym1;
end