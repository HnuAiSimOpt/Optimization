%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = L_bracketDomain(Demand,Arg)
  BdBox = [0 1 0 1];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),2/5);
  d2 = dRectangle(P,BdBox(1),2/5,BdBox(3),BdBox(4));
  Dist = dUnion(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  TopEdgeNodes = find(abs(Node(:,2)-BdBox(4))<eps);
  LowerNodes = find(abs(Node(:,2))<eps & ...
                    abs(Node(:,1))>0.9);
  RightNodes = find(abs(Node(:,1)-BdBox(2))<eps & ...
                    abs(Node(:,2))>0.75*2/5);
  RightUpperNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                    abs(Node(:,2)-2/5)<eps); 
  RightLowerNode = find(abs(Node(:,1)-BdBox(2))<eps & ...
                    abs(Node(:,2))<eps); 
  FixedNodes = TopEdgeNodes;
  Supp = zeros(length(FixedNodes),3);
  Supp(:,1) = FixedNodes; Supp(1:end,2) = 1; Supp(1:end,3) = 1;
  
  n1 = length(RightNodes); n2 = length(LowerNodes);
  n = n1+n2;
%   %Load = [RightNodes,zeros(n,1),-2*ones(n,1)/n]; % Load array 
  NStep = 200; % Number of time steps
  Th = linspace(0,pi,NStep+1);
  Load = zeros(n,2*(NStep+1)+1); 
  Load(1:n1,1) = RightNodes; Load(n1+1:n,1) = LowerNodes;
  Load(1:n1,3:2:end) = -sin(Th).*ones(n1,1)/n1; 
  Load(n1+1:n,2:2:end) = sin(Th).*ones(n2,1)/n2;
%   t = linspace(0,0.5,NStep+1);
%   Load = zeros(2,2*(NStep+1)+1);
%   Load(:,1) = [RightUpperNode;RightLowerNode];
% %   Load(1,3:2:end) = cos(4*pi*t+pi*1/3).*exp(-2*t)+1;
%   Load(1,3:2:end) = cos(Th).*exp(-2*t);
%   Load(2,3:2:end) = cos(Th).*exp(-2*t);
% % %   plot(Load(1,3:2:end));hold on
% %   plot(Load(2,3:2:end))
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%