function D = Pend3DModel_D(t,in2,in3)
%Pend3DModel_D
%    D = Pend3DModel_D(T,IN2,IN3)

%    This function was generated by the Symbolic Math Toolbox version 9.0.
%    08-Sep-2023 15:49:58

x1 = in2(:,1);
x2 = in2(:,2);
x3 = in2(:,3);
x4 = in2(:,4);
t2 = x2.^2;
t3 = x3.^2;
t4 = x4.^2;
t5 = x4.^3;
t6 = x1.*6.0e+1;
t7 = 1.0./x1;
t9 = x4.*1.2e+2;
t12 = x1.*x2.*x3.*1.2e+2;
t8 = -t6;
t10 = -t9;
t11 = t5.*1.2e+2;
t13 = t4.*x1.*1.2e+2;
t14 = t2.*t9;
t15 = t3.*t9;
D = reshape([(t7.*(t10+t11+t12+t14+t15))./1.0e+3,(t7.*(t8+t13+t3.*x1.*1.2e+2))./1.0e+2,0.0,0.0,0.0,0.0,t7.*(t8+t13+t2.*x1.*1.2e+2).*(-1.0./1.0e+3),(t7.*(t10+t11-t12+t14+t15))./1.0e+2,0.0,0.0,0.0,0.0,t7.*(x2.*-1.2e+2+t3.*x2.*1.2e+2+t4.*x2.*1.2e+2+x2.^3.*1.2e+2-x1.*x3.*x4.*1.2e+2).*(-1.0./1.0e+3),t7.*(x3.*-1.2e+2+t2.*x3.*1.2e+2+t4.*x3.*1.2e+2+x3.^3.*1.2e+2+t9.*x1.*x2).*(-1.0./1.0e+2),0.0,0.0,0.0,0.0],[6,3]);