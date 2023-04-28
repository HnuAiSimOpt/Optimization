%-------------------------- ConstraintsBridge ----------------------------%
VolFrac = 0.3;L = 30;
ElemCtrd = Centroids(fem);
Hb = max(Node(:,2)); % Height of bridge
fem.SElem = find(ElemCtrd(:,2)>=Hb-L/40); % Elements in passive region
ElemInd{1} = setdiff((1:fem.NElem)',fem.SElem);
%-------------------------------------------------------------------------%