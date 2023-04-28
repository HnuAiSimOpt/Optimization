%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = MbbDomain(Demand,Arg)
  BdBox = [0 3 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  Dist = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
  LeftUpperNode = find(abs(Node(:,1)-BdBox(1))<eps & ...
                       abs(Node(:,2)-BdBox(4))<eps);
  RightBottomNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                         abs(Node(:,2)-BdBox(3))<eps);
  RightUpperNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                         abs(Node(:,2)-BdBox(4))<eps);
  BottomMidNode = find(abs(Node(:,1)-0.5*BdBox(2))<eps & ...
                         abs(Node(:,2)-BdBox(3))<eps);
%   FixedNodes = [LeftEdgeNodes; RightBottomNode]; % Mbb  
%   Supp = zeros(length(FixedNodes),3);
%   Supp(:,1)=FixedNodes; Supp(1:end-1,2)=1; Supp(end,3)=1;   
%   Load = [LeftUpperNode,0,-0.5];
  FixedNodes = LeftEdgeNodes;
  Supp = ones(length(FixedNodes),3); Supp(:,1)=FixedNodes;
  NStep = 100; Th = linspace(0,pi,NStep+1);
  Load = zeros(2,2*(NStep+1)+1);
  Load(:,1) = [RightUpperNode;BottomMidNode];
  Load(1,3:2:end) = -sin(Th);
  Load(2,3:2:end) = cos(Th);
%   n = length(LeftUpperNode);
%   NStep = 200; % Number of time steps
%   Th = linspace(0,pi,NStep+1);
%   Load = zeros(n,2*(NStep+1)+1);
%   Load(:,1) = LeftUpperNode;
%   %Load(:,3:2:end) = -1e6*sin(Th).*ones(n,1)/n;
%   Load(:,3:2:end) = -1e6.*repmat(rand(1,NStep+1),n,1)/n;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%