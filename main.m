%% Define geometry
R = 1;
Mesh = Geo_Creation(R);
%% Integrate domain
[K,M] = Geo_Integration(Mesh);
%% Integrate source term
f =  @(x,y)-exp(-2*pi*x).*exp(-2*pi*y);
F = f(Mesh.XY(1,:)',Mesh.XY(2,:)');
F = -M*F;
%% Divide domain
Elms1 = find(all(ismember(Mesh.Topology,find((Mesh.XY(2,:)-Mesh.XY(1,:))>=0)),1));
Top1 = Mesh.Topology(:,Elms1);
Elms2 = setdiff(1:size(Mesh.Topology,2),Elms1);
Top2 = Mesh.Topology(:,Elms2);
% Plot subdomains
figure; patch('Faces',Top1(1:3,:)','Vertices',Mesh.XY','FaceColor','r','EdgeColor','k'); axis tight; axis equal;
patch('Faces',Top2(1:3,:)','Vertices',Mesh.XY','FaceColor','b','EdgeColor','k'); axis tight; axis equal;
%% Identify DOFs
N1 = unique(Top1);
N2 = unique(Top2);
Nb = intersect(N1,N2);
N1 = setdiff(N1,Nb);
N2 = setdiff(N2,Nb);
%% Constraint and Free DOFs
Nbf = setdiff(Nb,intersect(Nb,find(Mesh.Status)));
N1f = setdiff(N1,intersect(N1,find(Mesh.Status)));
N2f = setdiff(N2,intersect(N2,find(Mesh.Status)));
%% Isolate terms
K11 = K(N1f,N1f);
K22 = K(N2f,N2f);
K1b = K(N1f,Nbf);
Kb1 = K(Nbf,N1f);
K2b = K(N2f,Nbf);
Kb2 = K(Nbf,N2f);
Kbb = K(Nbf,Nbf);
F1 = F(N1f);
F2 = F(N2f);
Fb = F(Nbf);
%% Solve
u = zeros(numel(Mesh.Status),1);
%% Solve schur complement
A = Kbb-Kb1*(K11\K1b)-Kb2*(K22\K2b);
f = Fb-Kb1*(K11\F1)-Kb2*(K22\F2);
u(Nbf) = pcg(A,f,1e-9);
%% Solve domain 1
A = K11;
f = F1-K1b*u(Nbf);
u(N1f) = A\f;
%% Solve domain 2
A = K22;
f = F2-K2b*u(Nbf);
u(N2f) = A\f;
%% Plot result
figure; patch('Faces',Mesh.Topology(1:3,:)','Vertices',[Mesh.XY',u],'FaceVertexCData',u,'FaceColor','interp','EdgeColor','k'); axis tight; axis equal;
%% Solve reference
Rest_Nodes = find(Mesh.Status);
Free_Nodes = find(~Mesh.Status);
u = zeros(numel(Mesh.Status),1);
u(Free_Nodes) = K(Free_Nodes,Free_Nodes)\F(Free_Nodes);
figure; patch('Faces',Mesh.Topology(1:3,:)','Vertices',[Mesh.XY',u],'FaceVertexCData',u,'FaceColor','interp','EdgeColor','k'); axis tight; axis equal;