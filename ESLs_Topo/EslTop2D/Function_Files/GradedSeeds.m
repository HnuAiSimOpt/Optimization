%----------------------------- GradedSeeds -------------------------------%
function P = GradedSeeds(NElem,Domain,mu,Sym)
Pc = zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem  
  Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
  Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);                 %Index of seeds inside the domain
  prob = mu(Y(I,:)).^(1/2); prob = prob/max(prob);
  I = I(rand(size(I,1),1)<prob,:);
  NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
  Pc(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  Ctr = Ctr+NumAdded;
end
if ~exist('Sym','var')
  P = Pc; 
elseif strcmp(Sym,'X')==1
  P = [Pc(1:NElem/2,:);[Pc(1:NElem/2,1),-Pc(1:NElem/2,2)]]; 
elseif strcmp(Sym,'Y')==1
  P = [Pc(1:NElem/2,:);[-Pc(1:NElem/2,1),Pc(1:NElem/2,2)]]; 
elseif strcmp(Sym,'XY')==1
  P = [Pc(1:NElem/4,:);[-Pc(1:NElem/4,1),Pc(1:NElem/4,2)];...
      [-Pc(1:NElem/4,1),-Pc(1:NElem/4,2)];...
      [Pc(1:NElem/4,1),-Pc(1:NElem/4,2)]]; 
end
%-------------------------------------------------------------------------%