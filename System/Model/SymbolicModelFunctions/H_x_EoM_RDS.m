function H_x = H_x_EoM_RDS(in1)
%H_x_EoM_RDS
%    H_x = H_x_EoM_RDS(IN1)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    27-Jan-2024 19:58:07

u11 = in1(1,:);
u12 = in1(2,:);
u13 = in1(3,:);
t2 = u11.^2;
t3 = u12.^2;
t4 = u13.^2;
t5 = t2+t3+t4;
t6 = 1.0./t5;
t7 = t6.*u12.*u13.*2.745155490997198e+1;
t8 = -t7;
H_x = reshape([0.0,t6.*u11.*u12.*(-2.745155490997198e+1),t6.*u11.*u13.*(-2.745155490997198e+1),0.0,0.0,0.0,0.0,t3.*t6.*(-2.745155490997198e+1),t8,0.0,0.0,0.0,0.0,t8,t4.*t6.*(-2.745155490997198e+1),0.0,0.0,0.0],[6,3]);
