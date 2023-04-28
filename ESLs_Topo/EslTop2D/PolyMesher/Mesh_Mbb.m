%------------------------------ PolyDyna ---------------------------------%
% Ref: O Giraldo-Londoño, GH Paulino, "PolyDyna: A Matlab implementation  %
% for topology optimization of structures subjected to dynamic loads",    % 
% Structural and Multidisciplinary Optimization, 2021                     %
% DOI http://dx.doi.org/10.1007/s00158-021-02859-6                        %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load] = Mesh_Cantilever(Ne_ap)
L = 3; t = 1;
Ny = floor(round(sqrt(Ne_ap*t/L))); if mod(Ny,2)~=0, Ny = Ny+1; end
he = t/Ny; % Estimated element size (nn must be an even number)
Nx = Ny*L/t;
NElem = round(Nx*Ny); % Number of elements
[X,Y] = meshgrid(he/2:he:L-he/2,he/2:he:t-he/2); 
P = [X(:) Y(:)]; % Mesh seed
Domain = @MbbDomain;
[Node,Element,Supp,Load,~] = PolyMesher(Domain,NElem,0,P);
%-------------------------------------------------------------------------%