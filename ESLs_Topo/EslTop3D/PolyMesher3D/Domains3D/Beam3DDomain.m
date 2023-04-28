%-------------------------------------------------------------------------%
% PolyMesher3D developed by Anderson Pereira, April 2019                  %
% Contact: anderson@puc-rio.br                                            %
%-------------------------------------------------------------------------%
% Ref1: C Talischi, GH Paulino, A Pereira, IFM Menezes,                   %
%      "PolyMesher: A general-purpose mesh generator for polygonal        %
%      elements written in Matlab", Struct Multidisc Optim, 2012,         %
%      DOI 10.1007/s00158-011-0706-z                                      %
%                                                                         %
% Ref2: RS Thedin, A Pereira, IFM Menezes, GH Paulino,                    %
%      "Polyhedral mesh generation and optimization for finite element    %
%      computations. In: CILAMCE 2014 - XXXV Ibero-Latin American         %
%      Congress on Computational Methods in Engineering, 2014             %
%                                                                         %
% Ref3: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho,      %
%      "Implementation of fluid flow topology optimization in PolyTop",   %
%      Struct Multidisc Optim, 2016, DOI 10.1007/s00158-014-1182-z        %
%                                                                         %
% Ref4: H Chi, A Pereira, IFM Menezes, GH Paulino,                        %
%      "Virtual Element Method (VEM)-based topology optimization:         %
%      an integrated framework", Struct Multidisc Optim, 2019,            %
%      DOI 10.1007/s00158-019-02268-w                                     %
%-------------------------------------------------------------------------%
function [x] = Beam3DDomain(Demand,Arg)
  BdBox = [0 120 0 40 0 10];
  switch(Demand)
    case('Dist');   x = DistFnc(Arg,BdBox);
    case('BC');     x = BndryCnds(Arg{:},BdBox);
    case('BdBox');  x = BdBox;
    case('PFix');   x = FixedPoints(BdBox);
    case('Normal'); x = Normal(Arg{:});
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  d2 = [BdBox(5)-P(:,3), P(:,3)-BdBox(6)];
  d2 = [d2,max(d2,[],2)];
  Dist = dIntersect(d2,d1);
  %---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps =1e1*((BdBox(2)-BdBox(1))*...
             (BdBox(4)-BdBox(3))*...
             (BdBox(6)-BdBox(5)))^(1/3)/size(Element,1)/2;
  LeftNodes = find(abs(Node(:,1))<eps);
  RightNodes = find(abs(Node(:,1)-BdBox(2))<eps);
  SuppPart = [LeftNodes; RightNodes];
  Supp = ones(size(SuppPart,1),4);
  Supp(:,1) = SuppPart;
  UpperMiddleLineNodes = find(abs(Node(:,2)-BdBox(4))<eps & ...
            abs(Node(:,1)-0.5*BdBox(2))<2+eps);
  % Static Load
%   Load = -1*ones(size(LowerHalfCircleNodes,1),4);
%   Load(:,1) = LowerHalfCircleNodes; Load(:,2) = 0; Load(:,4) = 0;
  % Dynamic Load
  NStep = 100; % Number of time steps
  n = length(UpperMiddleLineNodes);
  Th = linspace(0,pi,NStep+1);
  Load = zeros(n,3*(NStep+1)+1);
  Load(:,1) = UpperMiddleLineNodes; Load(:,2) = 0; Load(:,4) = 0;
  Load(:,3:3:end) =-1e2*sin(Th).*ones(n,1)/2; % Half-cycle sinusoidal load
  
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
