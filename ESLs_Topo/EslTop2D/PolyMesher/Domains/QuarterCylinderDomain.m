%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = QuarterCylinder_Domain(Demand,Arg)
  BdBox = [0 0.2 0 0.2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  r1 = dRectangle(P,0,0.2,0,0.2);
  c1 = dCircle(P,BdBox(1),BdBox(3),BdBox(2));
  c2 = dCircle(P,BdBox(1),BdBox(3),BdBox(2)/2);
  R1 = dIntersect(c1,r1);
  Dist = dDiff(R1,c2);  
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  % Supports
  LeftEdgeNodes = find(abs(Node(:,1)-BdBox(1))<eps);
  BottomEdgeNodes = find(abs(Node(:,2)-BdBox(3))<eps);
  %FixedNodes = [LeftEdgeNodes; BottomEdgeNodes];
  FixedNodes = LeftEdgeNodes;
  Supp = zeros(length(FixedNodes),3); Supp(:,1)=FixedNodes;
  %Supp(1:length(LeftEdgeNodes),2)=1; Supp(1+length(LeftEdgeNodes):end,3)=1; 
  Supp(1:length(LeftEdgeNodes),2)=1;Supp(1:length(LeftEdgeNodes),3)=1;
  % Loads
  RightCornerNode = find(abs(Node(:,1)-(BdBox(2)))<eps & ...
                         abs(Node(:,2)-BdBox(3))<eps);
  InnerQuarterCircleNodes = find(sqrt(Node(:,1).^2+Node(:,2).^2)-0.1<eps);
  LowerNodes = find(Node(:,2)<eps & Node(:,1)-0.152<eps& Node(:,1)-0.148>eps); 
  n = length(LowerNodes);
  %Load = [RightCornerNode,400,0];
  NStep = 200; % Number of time steps
  Th = linspace(0,pi,NStep+1);
  Load = zeros(n,2*(NStep+1)+1); Load(:,1) = LowerNodes;
  Load(:,2:2:end) = -1e5*sin(Th).*ones(n,1)/n; 
  %Load(:,2:2:end) = -1e5.*rand(1,NStep+1).*ones(n,1)/n;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%