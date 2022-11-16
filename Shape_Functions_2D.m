function N = Shape_Functions_2D(X,Y,dFlag,EType,Degree)
if size(X,1)==1
    X = X';
end
if size(Y,1)==1
    Y = Y';
end
C0 = zeros(size(X));
C1 = ones(size(X));
switch EType
    case 1 % Triangle
        switch Degree
            case 1
                if dFlag
                    % Derivatives
                    dN1dXi = -1;
                    dN1dEta = -1;
                    dN2dXi = 1;
                    dN2dEta = 0;
                    dN3dXi = 0;
                    dN3dEta = 1;
                    N = [dN1dXi,dN2dXi,dN3dXi;dN1dEta,dN2dEta,dN3dEta];
                else
                    N1 = 1 - X - Y;
                    N2 = X;
                    N3 = Y;
                    N = [N1,N2,N3];
                end
            case 2
                if dFlag
                    dN1dXi = -3+4*X+4*Y;
                    dN1dEta = -3+4*X+4*Y;
                    dN2dXi = -1+4*X;
                    dN2dEta = C0;
                    dN3dXi = C0;
                    dN3dEta = -1+4*Y;
                    dN4dXi = +4-8*X-4*Y;
                    dN4dEta = -4*X;
                    dN5dXi = +4*Y;
                    dN5dEta = 4*X;
                    dN6dXi = -4*Y;
                    dN6dEta = +4-4*X-8*Y;
                    N = [dN1dXi,dN2dXi,dN3dXi,dN4dXi,dN5dXi,dN6dXi;dN1dEta,dN2dEta,dN3dEta,dN4dEta,dN5dEta,dN6dEta];
                else
                    N1 = (1-X-Y).*(1-2*X-2*Y);
                    N2 = -X.*(1-2*X);
                    N3 = -Y.*(1-2*Y);
                    N4 = 4*X.*(1-X-Y);
                    N5 = 4*X.*Y;
                    N6 = 4*Y.*(1-X-Y);
                    N = [N1,N2,N3,N4,N5,N6];
                end
        end
    case 2 % Square
        if dFlag
            % Derivatives
        else
        end
end
end