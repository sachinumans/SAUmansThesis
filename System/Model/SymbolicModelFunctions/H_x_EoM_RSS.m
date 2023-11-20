function H_x = H_x_EoM_RSS(in1)
%H_x_EoM_RSS
%    H_x = H_x_EoM_RSS(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    17-Nov-2023 12:43:05

u11 = in1(1,:);
u12 = in1(2,:);
u13 = in1(3,:);
t2 = u11.^2;
t3 = u12.^2;
t4 = u13.^2;
t5 = t2+t3+t4;
t6 = 1.0./t5;
t7 = t6.*u11.*u12.*2.745155490997198e+1;
t8 = t6.*u11.*u13.*2.745155490997198e+1;
t9 = t6.*u12.*u13.*2.745155490997198e+1;
t10 = -t7;
t11 = -t8;
t12 = -t9;
H_x = reshape([t2.*t6.*(-2.745155490997198e+1),t10,t11,0.0,0.0,0.0,t10,t3.*t6.*(-2.745155490997198e+1),t12,0.0,0.0,0.0,t11,t12,t4.*t6.*(-2.745155490997198e+1),0.0,0.0,0.0],[6,3]);
