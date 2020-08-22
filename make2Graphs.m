function [A,coor] = make2Graphs(N1, N2, type1, type2, C, varargin)

% Function that returns the adjacency matrix of the full system composed by 
% two graphs with "N1" and "N2" nodes and topologies "type1" and "type2", 
% respectively. Also returns the coordinates of its nodes. It also plots the
% resulting graph and computes the results using Schnakenberg model.

% INPUT:
%    [N1] = number of nodes of the first graph.
%    [N2] = number of nodes of the second graph.
%    [type1] and [typ2e] = string with the graphs topolgies: 1Dlattice, 
%                           2Dlattice, cycle, complete, regular, random.
%  [varargin] = are specified as parameter value pairs where the parameters are
%     'makeplot' = 1 -> plotting the graph topology
%     'circulant' = 1 x N vector with the circulant that generates the
%                   "regular" graph
%     'erdos' = scalar with the probability that generates the "random"
%               graph
%     'parameters' = 1 x 6 vector with the parameters of Schnakenberg model
%     'd' = 1 x 4 vector with the diffusion coefficients.
%     'edges'  = value for colouring the edges (1=black, 64=white)

% OUTPUT:
%  [A] = (N1+N2) x (N1+N2) matrix: adjacency matrix of the graph.
%  [coor] = (N1+N2) x 2 matrix: coordinates of the nodes for representation.


% Check graph topology 
if(length(type1)==length(type2))
    if(type1==type2 & (strcmp(type1,'regular') | (strcmp(type1,'random'))))
        [A1,coor1] = makeGraph(N1,type1,varargin{1:2});
        [A2,coor2] = makeGraph(N2,type2,varargin{3:4});
    else
        [A1,coor1] = makeGraph(N1,type1,varargin{:});
        [A2,coor2] = makeGraph(N2,type2,varargin{:});
    end
else
    [A1,coor1] = makeGraph(N1,type1,varargin{:});
    [A2,coor2] = makeGraph(N2,type2,varargin{:});
end

% Read values of the optional parameters:
nVarArgin=length(varargin);
for k=1:2:nVarArgin
    switch lower(varargin{k})
        case 'parameters'
            param = varargin{k+1};
        case 'd'
            d = varargin{k+1};
        case 'e'
            e = varargin{k+1};
        case 'edges'
            edg = varargin{k+1};
    end
end

% Construct global adjacency matrix and coordinates.
A = zeros(N1+N2);
A(1:N1,1:N1) = A1; A(N1+1:end, N1+1:end) = A2; A(1:N1,N1+1:end) = e*C; A(N1+1:end,1:N1) = e*C';  
coor1(:,1) = coor1(:,1);
coor2(:,2) = coor2(:,2)+2;
coor = [coor1;coor2];
wgPlot2(A, coor,'vertexMetadata',ones(N1+N2,1),'vertexScale',200,'edgeWidth',2,'edges',1);


% If the parameters are provided we solve the system of ODEs and plot the
% results.
if exist('param','var')
    [U1, U2, V1, V2, evalFull] = solve2Graphs(A1, A2, C, param, d, e);
    wgPlot2(A, coor,'vertexMetadata',[U1(end,:),U2(end,:)],'vertexScale',200,'edgeWidth',2,'edges',edg);
    caxis([0.6, 1]);
end


end