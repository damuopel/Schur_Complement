function [K,M] = Geo_Integration(Mesh)
% Constants
DOFsPerNode = 1;
NodesPerElement = size(Mesh.Topology,1);
% Integration quadrature
[xGauss,wGauss] = Quadrature_Library(2,1);
% Mesh
Topology = Mesh.Topology;
XY = Mesh.XY;
nElms = size(Topology,2);
% Initialize K 
K = zeros((DOFsPerNode*NodesPerElement)^2,nElms);
% Initialize M
M = zeros((DOFsPerNode*NodesPerElement)^2,nElms);
for iElm = 1:nElms
    iCoords = XY(:,Topology(:,iElm));
    for iGP = 1:numel(wGauss)
        Xi = xGauss(iGP,1);
        Eta = xGauss(iGP,2);
        % K Matrix
        dNlocal = Shape_Functions_2D(Xi,Eta,1,1,NodesPerElement/3);
        Jacobian = dNlocal*iCoords';
        dNGlobal = Jacobian\dNlocal;
        Ke = dNGlobal'*dNGlobal*wGauss(iGP)*det(Jacobian);
        K(:,iElm) = K(:,iElm) + Ke(:);
        % M Matrix
        N = Shape_Functions_2D(Xi,Eta,0,1,NodesPerElement/3);
        Me = N'*N*wGauss(iGP)*det(Jacobian);
        M(:,iElm) = M(:,iElm) + Me(:);
    end
end
iK = kron(Topology,ones(1,DOFsPerNode*NodesPerElement));
jK = kron(Topology,ones(DOFsPerNode*NodesPerElement,1));
% Assemble K Matrix
K = sparse(iK(:),jK(:),K(:));
% Assemble M Matrix
M = sparse(iK(:),jK(:),M(:));
end