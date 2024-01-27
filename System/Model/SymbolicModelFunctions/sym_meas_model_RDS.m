function y = sym_meas_model_RDS(in1,in2,in3)
%sym_meas_model_RDS
%    Y = sym_meas_model_RDS(IN1,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    27-Jan-2024 19:58:08

bS1 = in3(1,:);
bS2 = in3(2,:);
bS3 = in3(3,:);
u11 = in2(1,:);
u12 = in2(2,:);
u13 = in2(3,:);
u21 = in2(4,:);
u22 = in2(5,:);
u23 = in2(6,:);
u24 = in2(7,:);
u31 = in2(8,:);
u32 = in2(9,:);
u33 = in2(10,:);
u34 = in2(11,:);
u41 = in2(12,:);
u42 = in2(13,:);
u43 = in2(14,:);
u44 = in2(15,:);
x1 = in1(1,:);
x2 = in1(2,:);
x3 = in1(3,:);
t2 = u11.*x1;
t3 = u12.*x2;
t4 = u13.*x3;
t5 = u11.^2;
t6 = u12.^2;
t7 = u13.^2;
t8 = u21.*u32.*2.0;
t9 = u22.*u31.*2.0;
t10 = u21.*u33.*2.0;
t11 = u23.*u31.*2.0;
t12 = u21.*u34.*2.0;
t13 = u22.*u33.*2.0;
t14 = u23.*u32.*2.0;
t15 = u24.*u31.*2.0;
t16 = u22.*u34.*2.0;
t17 = u24.*u32.*2.0;
t18 = u23.*u34.*2.0;
t19 = u24.*u33.*2.0;
t20 = u21.*u42.*2.0;
t21 = u22.*u41.*2.0;
t22 = u21.*u43.*2.0;
t23 = u23.*u41.*2.0;
t24 = u21.*u44.*2.0;
t25 = u22.*u43.*2.0;
t26 = u23.*u42.*2.0;
t27 = u24.*u41.*2.0;
t28 = u22.*u44.*2.0;
t29 = u24.*u42.*2.0;
t30 = u23.*u44.*2.0;
t31 = u24.*u43.*2.0;
t32 = -t9;
t33 = -t11;
t34 = -t13;
t35 = -t15;
t36 = -t17;
t37 = -t18;
t38 = -t21;
t39 = -t23;
t40 = -t25;
t41 = -t27;
t42 = -t29;
t43 = -t30;
t44 = t2+t3+t4;
t45 = t5+t6+t7;
t46 = sqrt(t45);
t48 = t8+t19+t32+t37;
t49 = t10+t16+t33+t36;
t50 = t12+t14+t34+t35;
t51 = t20+t31+t38+t43;
t52 = t22+t28+t39+t42;
t53 = t24+t26+t40+t41;
t47 = 1.0./t46;
t54 = bS1.*t49;
t55 = bS1.*t50;
t56 = bS2.*t48;
t57 = bS2.*t50;
t58 = bS3.*t48;
t59 = bS3.*t49;
t63 = t46.*7.07125314117414e+3;
t60 = -t56;
t61 = -t58;
t62 = -t59;
t64 = -t63;
t68 = t44.*t47.*1.57297407539752e+3;
t65 = t54+t60;
t66 = t55+t61;
t67 = t57+t62;
t69 = t64+t68+7.311756864837817e+3;
mt1 = [-bS2.*t53+bS3.*t52-t49.*t65-t50.*t66+u21.*u23.*3.924e+1-u22.*u24.*3.924e+1];
mt2 = [bS1.*t53-bS3.*t51+t48.*t65-t50.*t67-u21.*u22.*3.924e+1-u23.*u24.*3.924e+1-t47.*t69.*u12.*1.745200721317321e-2];
mt3 = [-bS1.*t52+bS2.*t51+t48.*t66+t49.*t67-u21.^2.*1.962e+1+u22.^2.*1.962e+1+u23.^2.*1.962e+1-u24.^2.*1.962e+1-t47.*t69.*u13.*1.745200721317321e-2;t48;t49;t50];
y = [mt1;mt2;mt3];
