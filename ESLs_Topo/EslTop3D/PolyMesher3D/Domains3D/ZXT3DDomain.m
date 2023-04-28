%-------------------------------------------------------------------------%
% 轴箱体                                   %
%-------------------------------------------------------------------------%
function [x] = ZXT3DDomain(Demand,Arg)
  BdBox = [-90, 720, -170, 160, -115, 115];
  switch(Demand)
    case('Dist');   x = DistFnc(Arg,BdBox);
    case('BC');     x = BndryCnds(Arg{:},BdBox);
    case('BdBox');  x = BdBox;
    case('PFix');   x = FixedPoints(BdBox);
    case('Normal'); x = Normal(Arg{:});
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  c1 = dCircle(P,0,0,70);
  c2 = dCircle(P,0,0,90);
  c3 = dCircle(P,450,0,120);
  c4 = dCircle(P,450,0,140);
%   N_p = [45,45;45,-45;125,-25;235,-25;270,-45;270,-65;
%          310,-65;350,-110;350,-150;550,-150;550,-110;
%          570,-80;610,-80;610,-45;720,-45;720,0;610,0;
%          610,30;565,75;565,115;335,115;155,45];
  N_p = [45,77.9423;45,-77.9423;350,-170;550,-170;
         720,-45;720,0;565,160;335,160];
  d1 = dPolygon(P,N_p);
  douter = dUnion(c2,d1);
  din = dUnion(c1,c3);
  Dist = dDiff(douter,din);
  
  N_p1 = [610,0;720,0;565,160];
  d3 = dPolygon(P,N_p1);
  Dist = dDiff(Dist,d3);
  % Extrusion
  d2 = [BdBox(5)-P(:,3), P(:,3)-BdBox(6)];
  d2 = [d2,max(d2,[],2)];
  Dist = dIntersect(d2,Dist);
  %---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps =1e1*((BdBox(2)-BdBox(1))*...
             (BdBox(4)-BdBox(3))*...
             (BdBox(6)-BdBox(5)))^(1/3)/size(Element,1)/2;
  RightCircleNodes = find(abs(sqrt((Node(:,1)-450).^2+Node(:,2).^2))<120+eps);
  Supp = ones(size(RightCircleNodes,1),4);
  Supp(:,1) = RightCircleNodes;
  LeftCircleNodes = ...
      find(abs(sqrt(Node(:,1).^2+Node(:,2).^2))<70+eps);
  UpperNodes = find(abs(Node(:,2)-BdBox(4))<eps);
  UpperAnnulusNodes = find(abs(Node(:,2)-BdBox(4))<eps & ...
      abs(sqrt((Node(:,1)-450).^2+Node(:,3).^2))<100+eps & ...
      abs(sqrt((Node(:,1)-450).^2+Node(:,3).^2))>50+eps);
  RightHoleNode1 = ...
      sqrt((Node(:,1)-660).^2+Node(:,2).^2+(Node(:,3)+80).^2);
  [foo,RightHoleNode1] = sort(RightHoleNode1);
  RightHoleNode2 = ...
      sqrt((Node(:,1)-660).^2+Node(:,2).^2+(Node(:,3)-20).^2);
  [foo,RightHoleNode2] = sort(RightHoleNode2);
  RightHoleNodes = [RightHoleNode1(1);RightHoleNode2(1)];
  % Dynamic Load
  NStep = 100; % Number of time steps
  n1 = length(LeftCircleNodes);
  n2 = length(UpperAnnulusNodes);
  n3 = length(RightHoleNodes);
  n = n1+n2+n3;
  Th = linspace(0,pi,NStep+1);
  Load = zeros(n,3*(NStep+1)+1);
  
  Load(1:n1,1) = LeftCircleNodes;
  Load(n1+1:n1+n2,1) = UpperAnnulusNodes;
  Load(n1+n2+1:n,1) = RightHoleNodes;
  
  Load(1:n1,2:3:end) = 1e2*sin(Th).*ones(n1,1)/n1; 
  Load(1:n1,4:3:end) = 1e2*sin(Th).*ones(n1,1)/n1; 
  
  Load(n1+1:n1+n2,3:3:end) = -1e2*sin(Th).*ones(n2,1)/n2;
  
  Load(n1+n2+1:n,3:3:end) = -1e2*sin(Th).*ones(n3,1)/n3;
%   % Static Load
%   Load = zeros(n,4);
%   Load(1:n1,1) = LeftCircleNodes;
%   Load(n1+1:n1+n2,1) = UpperAnnulusNodes;
%   Load(n1+n2+1:n,1) = RightHoleNodes;
%   
%   Load(1:n1,2) = 1e2.*ones(n1,1)/n1; 
%   Load(1:n1,4) = 1e2.*ones(n1,1)/n1; 
%   
%   Load(n1+1:n1+n2,3) = -1e2.*ones(n2,1)/n2;
%   
%   Load(n1+n2+1:n,3) = -1e2.*ones(n3,1)/n3;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%--------------------------------------------- SPECIFY ANALYTICAL GRADIENTS
function [x] = Normal(P,n1,n2,n3)
  n1(:,1:2) = 0;
  n2(:,1:2) = 0;
  n3(:,1) = -1; n3(:,2) =  1;
  n3(:,3:end) =  0;
  x = {n1,n2,n3};
%-------------------------------------------------------------------------%
