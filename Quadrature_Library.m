function [Coords,Weights] = Quadrature_Library(Dimension,EType)
% Dimension -> 1 (1D), 2 (2D) and 3 (3D)
% EType -> Element type: 1 (Triangle/Tetrahedron) and 2 (Square/Hexahedron)
switch Dimension
    case 1
        Coords=[...
            -0.577350269189626;...
            0.577350269189626];
        Weights=[...
            1;...
            1];
    case 2
        switch EType
            case 1 % Triangle
                Coords=[...
                    1/6,1/6;...
                    2/3,1/6;...
                    1/6,2/3];
                Weights=[...
                    1/6;...
                    1/6;...
                    1/6];
            case 2 % Square
                Coords=[...
                    -0.57735026918963, -0.57735026918963;...
                    +0.57735026918963, -0.57735026918963;...
                    -0.57735026918963, +0.57735026918963;...
                    +0.57735026918963  +0.57735026918963];
                Weights=[...
                    1;...
                    1;...
                    1;...
                    1];
        end
    case 3
        switch EType
            case 1 % Tetrahedron
                Coords=[...
                    -0.7236067977 -0.7236067977 -0.7236067977
                    0.1708203932 -0.7236067977 -0.7236067977
                    -0.7236067977 0.1708203932 -0.7236067977
                    -0.7236067977 -0.7236067977 0.1708203932 ];
                Weights =[1/3 1/3 1/3 1/3 ];
            case 2 % Hexahedron
                Coords=[...
                    -0.57735026918963, -0.57735026918963, -0.57735026918963;...
                    +0.57735026918963, -0.57735026918963, -0.57735026918963;...
                    -0.57735026918963, +0.57735026918963, -0.57735026918963;...
                    +0.57735026918963  +0.57735026918963, -0.57735026918963;...
                    -0.57735026918963, -0.57735026918963, +0.57735026918963;...
                    +0.57735026918963, -0.57735026918963, +0.57735026918963;...
                    -0.57735026918963, +0.57735026918963, +0.57735026918963;...
                    +0.57735026918963  +0.57735026918963, +0.57735026918963];
                Weights=[...
                    1;...
                    1;...
                    1;...
                    1;...
                    1;...
                    1;...
                    1;...
                    1];
        end
end
end