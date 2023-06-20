
load ExplEoM_LSS.mat 
vSS = {[x1 x2 x3 x4 x5 x6 x7 x8 x9 x10 x11 x12 x13 x14], t, [F1,F2], Vl,Vs,h,Wi,l0,m, k,b, [j_in1,j_in2,j_in3]};
matlabFunction(dx,'File','LSSeom','Vars', vSS)

load ExplEoM_RSS.mat 
matlabFunction(dx,'File','RSSeom','Vars', vSS)

load ExplEoM_DS.mat
matlabFunction(dx,'File','DSeom')
