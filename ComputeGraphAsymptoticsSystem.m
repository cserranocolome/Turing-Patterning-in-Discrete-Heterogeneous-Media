function [evalFull, evecFull, evalAsym0, evecAsym0, evalAsym1, lambdaAsym] = ComputeGraphAsymptoticsSystem(A1, A2, C0, e, J1, J2, d)

% Function that returns the asymptotic approximations of the eigenvalues of
% the full system composed by two coupled graphs. System case.

% INPUT:
%    [A1] = N1 x N1 matrix. Adjacency matrix of the first graph.
%    [A2] = N2 x N2 matrix. Adjacency matrix of the second graph.
%  [C0] = N1 x N2 matrix. Coupling matrix
%   [e] = epsilon (strength of the coupling)
%   [J1] = 2 x 2 Jacobian matrix of the first graph
%   [J2] = 2 x 2 Jacobian matrix of the second graph
%   [d] = 1 x 4 vector containin the diffusion coefficients: du1,dv1,du2,dv2
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

du1 = d(1);
dv1 = d(2);
du2 = d(3);
dv2 = d(4);

fu1 = J1(1,1); fv1 = J1(1,2); gu1 = J1(2,1); gv1 = J1(2,2);
fu2 = J2(1,1); fv2 = J2(1,2); gu2 = J2(2,1); gv2 = J2(2,2);

% Create two lattices
%L1 = A1-diag(sum(A1,1));
%L2 = A2-diag(sum(A2,1));

D1 = -diag(sum(C0,2));
D2 = -diag(sum(C0,1));

M1 = zeros(2*N1+2*N2);  M1(1:N1,1:N1) = du1*D1; M1(1:N1,2*N1+1:2*N1+N2) = du1*C0; M1(N1+1:2*N1,N1+1:2*N1) = dv1*D1; M1(N1+1:2*N1,2*N1+N2+1:end) = dv1*C0;
M1(2*N1+1:2*N1+N2,1:N1) = du2*C0'; M1(2*N1+1:2*N1+N2,2*N1+1:2*N1+N2) = du2*D2; M1(2*N1+N2+1:end,N1+1:2*N1) = dv2*C0'; M1(2*N1+N2+1:end,2*N1+N2+1:end) = dv2*D2;


%Full Graph
C = e*C0;
A0 = [A1,0*C;0*C',A2];
A3 = [A1,C;C',A2];

%Full Laplacian
L0 = A0-diag(sum(A0,1));
L1 = A1-diag(sum(A1,1));
L2 = A2-diag(sum(A2,1));
L3 = A3-diag(sum(A3,1));

%L3-L0-(D*e+[0*A1,C;C',0*A2])


%Full system
M11 = fu1*eye(N1) + du1*L1; M12 = fv1*eye(N1); M21 = gu1*eye(N1); M22 = gv1*eye(N1)+dv1*L1;
M33 = fu2*eye(N2) + du2*L2; M34 = fv2*eye(N2); M43 = gu2*eye(N2); M44 = gv2*eye(N2)+dv2*L2;
M0 = zeros(2*N1+2*N2); M0(1:N1,1:N1) = M11; M0(1:N1,N1+1:2*N1) = M12; M0(N1+1:2*N1,1:N1) = M21; M0(N1+1:2*N1,N1+1:2*N1) = M22;
M0(2*N1+1:2*N1+N2,2*N1+1:2*N1+N2) = M33; M0(2*N1+1:2*N1+N2,2*N1+N2+1:end) = M34; M0(2*N1+N2+1:end,2*N1+1:2*N1+N2) = M43; M0(2*N1+N2+1:end,2*N1+N2+1:end) = M44;

M = M0+M1*e;    %J+L+extra connections
%M0 is J+L

%Eigenvalues and eigenvectors of the full system
[VFull,evalFull] = eig(M);
[evalFull,indFull] = esort(diag(evalFull));
evecFull = VFull(:,indFull);

% Asymptotic evalues:
[V,evalAsym0,W] = eig(M0);
%[V2,evalAsym02] = eig(M0');
[evalAsym0,ind] = esort(diag(evalAsym0));
%[evalAsym02,ind2] = esort(diag(evalAsym02));
evecAsym0 = V(:,ind);
W = W(:,ind);
%evalAsym1 = diag(evecAsym0'*M1*evecAsym0);
%evalAsym1 = diag(conj(V2(:,ind2))'*M1*evecAsym0)./diag(conj(V2(:,ind2))'*evecAsym0);
evalAsym1 = diag(W'*M1*evecAsym0)./diag(W'*evecAsym0);
%conj(V2(:,ind2))'*evecAsym0

lambdaAsym = evalAsym0 + e*evalAsym1;
end