%----- ------------------------- PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = EyeBarDomain(Demand,Arg)
%   BdBox = [0 1.6 -0.4 0.4];
  BdBox = [0 2 0 1];
%   BdBox = [0 2000 0 1000];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
%   d2 = dCircle(P,3*BdBox(4),0,0.15);
  d2 = dCircle(P,1,0.5,0.3);
%   d2 = dCircle(P,1000,500,300);
  Dist = dDiff(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  R = 0.15; C = BdBox(4);
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
  BottomMidNode = sqrt((Node(:,1)-0.5*BdBox(2)).^2+(Node(:,2)-BdBox(3)).^2);
  [foo,BottomMidNode] = sort(BottomMidNode);
  RightTopNode = sqrt((Node(:,1)-BdBox(2)).^2+(Node(:,2)-BdBox(4)).^2);
  [foo,RightTopNode] = sort(RightTopNode);
  Supp = ones(length(LeftEdgeNodes),3); Supp(:,1) = LeftEdgeNodes;
  NStep = 200; % Number of time steps
  Th = linspace(0,pi,NStep+1);
  Load = zeros(2,2*(NStep+1)+1); 
  Load(:,1) = [RightTopNode(1);BottomMidNode(1)];
  Load(1,3:2:end) = -sin(Th); Load(2,3:2:end) = cos(Th);
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%