nelx = 60; 
nely = 20; 
nelz = 12; 
volfrac = 0.3; 
er = 0.04; 
rmin = 1.5; 
start = 25;
[TFE,Comp,err,loop,Tiji] = MgcgRomETO3d(nelx,nely,nelz,volfrac,er,rmin,start);