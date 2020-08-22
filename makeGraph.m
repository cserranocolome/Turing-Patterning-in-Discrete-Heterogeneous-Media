function [A,coor] = makeGraph(N,type,varargin)

% Function that returns the adjacency matrix of a graph with "N" nodes and
% topology "type" and the coordinates of its nodes. It also plots the
% resulting graph and computes the results using Schnakenberg model.

% INPUT:
%    [N] = number of nodes 
%     [type] = string with the graph topolgy: 1Dlattice, 2Dlattice, cycle, complete, regular, random.
%  [varargin] = are specified as parameter value pairs where the parameters are
%     'makeplot' = 1 -> plotting the graph topology
%     'model' = 1 -> solves the system and plots U
%     'circulant' = 1 x N vector with the circulant that generates the
%                   "regular" graph
%     'erdos' = scalar with the probability that generates the "random"
%               graph
%     'parameters' = 1 x 3 vector with the parameters of Schnakenberg model
%     'd' = 1 x 2 vector with the diffusion coefficients.
%       'edges'  = value for colouring the edges (1=black, 64=white)

% OUTPUT:
%  [A] = N x N matrix: adjacency matrix of the graph.
%  [coor] = N x 2 matrix: coordinates of the nodes for representation.

% Read values of the optional parameters:
makeplot = 0;
nVarArgin=length(varargin);
for k=1:2:nVarArgin
    switch lower(varargin{k})
        case 'circulant'
            circ = varargin{k+1};
        case 'erdos'
            p = varargin{k+1};
        case 'makeplot'
            makeplot = 1;
        case 'parameters'
            param = varargin{k+1};
        case 'd'
            d = varargin{k+1};
        case 'model'
            model = varargin{k+1};
        case 'edges'
            edg = varargin{k+1};
    end
end

% Check the graph topology
switch type
    case '1Dlattice'
        A = diag(ones(N-1,1),1); A = A+A';
        coor = graphplot(A,'1Dlattice');
        
        if makeplot==1
            wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
        end
    case '2Dlattice'
        if mod(sqrt(N),1)~=0
            error('N needs to be a square number for a 2D lattice');
        else
            I = eye(sqrt(N));
            A1 = diag(ones(sqrt(N)-1,1),1); A1 = A1+A1';
            A = kron(I,A1)+diag(ones(sqrt(N)*(sqrt(N)-1),1),sqrt(N))+diag(ones((sqrt(N)-1)*sqrt(N),1),-sqrt(N));
            coor = graphplot(A,'2Dlattice');
            if makeplot==1
                wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
            end
        end
    case 'cycle'
        A = diag(ones(N-1,1),1); A(1,N) = 1; A = A+A';
        coor = graphplot(A,'cycle');
        if makeplot==1
            wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
        end
    case 'complete'
        A = (triu(ones(N))-eye(N)); A = A+A';
        coor = graphplot(A,'cycle');
        if makeplot==1
            wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
        end
    case 'star'
        A = zeros(N); A(1,2:end) = ones(N-1,1); A = A+A';
        coor = graphplot(A,'wheel');
        if makeplot==1
            wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
        end
    case 'regular'
        A = circulant(circ,1); 
        coor = graphplot(A,'cycle');
        if makeplot==1
            wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
        end
    case 'random'
        A = erdos_reyni(N,p); A = make_connected(A); 
        coor = graphplot(A,'cycle');
        if makeplot==1
            wgPlot2(A, coor,'vertexMetadata',ones(N,1),'vertexScale',200,'edgeWidth',2,'black','k');
        end
    otherwise
        fprintf('Error, no such type is found!\n');
end

% If the model parameter is provided we solve the system of ODEs and plot the
% results.

if exist('model','var')
    [U1, V1] = solve1Graph(A, param, d);
    wgPlot2(A, coor,'vertexMetadata',U1(end,:),'vertexScale',200,'edgeWidth',2);
    caxis([0.1, 1.1]);
end

end