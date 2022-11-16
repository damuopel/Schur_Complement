function Mesh = Geo_Creation(R)
TOL_GEO = 1e-2;
Num_Radial_Elms = 10;
Num_Angle_Elms = 30;
r = linspace(0,R,Num_Radial_Elms);
t = linspace(0,2*pi,Num_Angle_Elms);
% Point Cloud
C = [0;0];
X = C(1)+r'*cos(t);
Y = C(2)+r'*sin(t);
P = [X(:),Y(:)];
DT = delaunayTriangulation(P);
% Create Mesh
model = createpde;
geometryFromMesh(model,DT.Points',DT.ConnectivityList');
auxMesh = generateMesh(model,'GeometricOrder','linear','Hmax',R/10);
% auxMesh = generateMesh(model,'GeometricOrder','quadratic','Hmax',R/10);
Mesh = struct;
Mesh.XY = auxMesh.Nodes;
Mesh.Topology = auxMesh.Elements;
Mesh.Status = abs(sqrt(sum((Mesh.XY-C).^2,1))-R)<=TOL_GEO;
% scatter(Mesh.XY(1,:),Mesh.XY(2,:),20,Mesh.Status);
% patch('Faces',Mesh.Topology(1:3,:)','Vertices',Mesh.XY','FaceColor','r','EdgeColor','k');
end