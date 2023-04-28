%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = HornDomain(Demand,Arg)
  BdBox = [-1 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dCircle(P,0,0,1);
  d2 = dCircle(P,-0.4,0,0.55);
  d3 = dLine(P,0,0,1,0);
  Dist = dIntersect(d3,dDiff(d1,d2));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  LowerEdgeNodes = find(abs(Node(:,2))<eps&Node(:,1)<0);
  LowerEdgeNodes1 = find(abs(Node(:,2))<eps&Node(:,1)>0);
  RightNode = find(abs(Node(:,1)-BdBox(2))<eps&abs(Node(:,2))<eps);
  TopNodes = sqrt(Node(:,1).^2+(Node(:,2)-BdBox(4)).^2);
  [foo,TopNodes] = sort(TopNodes);
  TopNode = TopNodes(1);
  FixedNodes = LowerEdgeNodes;
  Supp = ones(length(FixedNodes),3);
  Supp(:,1)=FixedNodes; 
  
  n = length(LowerEdgeNodes1);
  NStep = 200; % Number of time steps
  Th = linspace(0,pi,NStep+1);
  Load = zeros(n,2*(NStep+1)+1);
  Load(:,1) = LowerEdgeNodes1;
  Load(:,3:2:end) = 1e6*sin(Th).*ones(n,1)/n;
  %Load(:,3:2:end) = -1e6.*repmat(rand(1,NStep+1),n,1)/n;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%