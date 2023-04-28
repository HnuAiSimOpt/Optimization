%-------------------------------- Areas ----------------------------------%
function ElemArea = Areas(fem)
ElemArea = zeros(fem.NElem,1); 
for el = 1:fem.NElem 
  vx = fem.Node(fem.Element{el},1); vy = fem.Node(fem.Element{el},2);
  temp = vx.*vy([2:end 1])-vy.*vx([2:end 1]);
  ElemArea(el,1) = 0.5*sum(temp);
end
%-------------------------------------------------------------------------%