%------------------------- ConstraintsBuilding ---------------------------%
NConstr = 5; 
B = 30; H = 75;
Hc = H/NConstr;
ElemCtrd = Centroids(fem);
ElemInd = cell(NConstr,1);
for ii=1:NConstr
  ElemInd{ii} = find(abs(ElemCtrd(:,2)-(ii-1/2)*Hc)<=Hc/2);
end
sc = 0.8; a = sc*2/3*H; b = sc*B/2; c = sc*(H-a);
fem.SElem = find((((ElemCtrd(:,2)-c)/a).^2+(ElemCtrd(:,1)/b).^2>=...
                1+30*eps)&(ElemCtrd(:,2)<=4/5*H));
for ii=1:NConstr; ElemInd{ii} = setdiff(ElemInd{ii},fem.SElem); end
ElemArea = Areas(fem);
Vmax = 0.5*sum(ElemArea(ElemInd{5}));
for ii=1:NConstr; VolFrac(ii,1) = Vmax/sum(ElemArea(ElemInd{ii})); end
% VolFrac(end,1) = 0.75;
%-------------------------------------------------------------------------%