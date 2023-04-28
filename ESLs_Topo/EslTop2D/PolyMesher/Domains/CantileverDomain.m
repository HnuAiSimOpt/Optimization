%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
%% RectangleDomain
function [x] = CantileverDomain(Demand,Arg)
  BdBox = [0 8 -2 2];
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
  RightEdgeNodes = find(abs(Node(:,1)-BdBox(2))<eps);
  RigthCenterNode = find(abs(Node(:,1)-BdBox(2))<eps&abs(Node(:,2))<eps);
  LeftCenterNode = find(abs(Node(:,1))<eps&abs(Node(:,2))<eps);
  RigthBottomNode = find(abs(Node(:,1)-BdBox(2))<eps&abs(Node(:,2)-BdBox(3))<eps);
  UpperCenterNode = find(abs(Node(:,1)-1/2*BdBox(2))<eps&abs(Node(:,2)-BdBox(4))<eps);
  % Another
%   UpperCenterNode = sqrt((Node(:,1)-(BdBox(2)/2)).^2+...
%                          (Node(:,2)-BdBox(4)).^2);
%   [foo,UpperCenterNode] = sort(UpperCenterNode);
  
  InnerNode1 = find(abs(Node(:,1)-1/4*BdBox(2))<eps&abs(Node(:,2))<eps);
  InnerNode2 = find(abs(Node(:,1)-3/4*BdBox(2))<eps&abs(Node(:,2))<eps);
  
  FixedNodes = LeftEdgeNodes;
  Supp = ones(length(FixedNodes),3);  Supp(:,1)=FixedNodes; 
  % case1:CantileverBeam
  NStep = 100; % Number of time steps
  Th = linspace(0,pi,NStep+1);
  Load = zeros(1,2*(NStep+1)+1); Load(1,1) = RigthCenterNode;
  Load(1,3:2:end) = -1000*sin(Th); % Half-cycle sinusoidal load
  % case2:Supp = [RigthCenterNode;LeftCenterNode]
%   n1 = length(InnerNode1); n2 = length(InnerNode2);
%   n = n1+n2;
%   NStep = 200; % Number of time steps
%   Th = linspace(0,pi,NStep+1);
%   Load = zeros(n,2*(NStep+1)+1);
%   Load(1:n1,1) = InnerNode1;
%   Load(n1+1:end,1) = InnerNode2;
%   Load(1:n1,3:2:end) = 1e6*sin(Th).*ones(n1,1)/n1;
%   Load(n1+1:end,3:2:end) = -1e6*sin(Th).*ones(n2,1)/n2;
  % case3:Supp = [LeftEdgeNodes;RightEdgeNodes]  
%   n = length(UpperCenterNode);
%   NStep = 200; % Number of time steps
%   Th = linspace(0,pi,NStep+1);
%   Load = zeros(n,2*(NStep+1)+1);
%   Load(:,1) = UpperCenterNode;
%   Load(:,3:2:end) = -1e6*sin(Th).*ones(n,1)/n;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%