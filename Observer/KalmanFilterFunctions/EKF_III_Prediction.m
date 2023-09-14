function [m_mink,P_mink] = EKF_III_Prediction(f, m_km, u_km, F_x, F_xx, P_km, Q)
%EKF_I_PREDICTION Extended Kalman filter variant 3
%   Sarkka book page 73
nx = length(m_km);
SUMm = zeros(nx, 1);
SUMp = zeros(nx);
for i = 1:nx
    e = [zeros(i-1, 1); 1; zeros(nx-i, 1)];
    SUMm = SUMm + e*trace(F_xx(:,:,i)*P_km);
    for i2 = 1:nx
        e2 = [zeros(i2-1, 1); 1; zeros(nx-i2, 1)];
        SUMp = SUMp + e*e2'*trace(F_xx(:,:,i)*P_km*F_xx(:,:,i2)*P_km);
    end
end

m_mink = f(0, m_km', u_km) + 0.5*SUMm;
P_mink = F_x*P_km*F_x' + 0.5*SUMp + Q;
end

