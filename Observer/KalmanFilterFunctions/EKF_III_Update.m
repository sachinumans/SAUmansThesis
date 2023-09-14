function [m_k, P_k] = EKF_III_Update(y, h, m_mink, u_k, H_x, H_xx, P_mink, R, frobNormBound, rescaleFactor)
%EKF_I_UPDATE Extended Kalman filter variant 3
%   Sarkka book page 73
ny = length(y);
SUMv = zeros(ny, 1);
SUMs = zeros(ny);
for i = 1:ny
    e = [zeros(i-1, 1); 1; zeros(ny-i, 1)];
    SUMv = SUMv + e*trace(H_xx(:,:,i)*P_mink);
    for i2 = 1:ny
        e2 = [zeros(i2-1, 1); 1; zeros(ny-i2, 1)];
        SUMs = SUMs + e*e2'*trace(H_xx(:,:,i)*P_mink*H_xx(:,:,i2)*P_mink);
    end
end


v = y - h(0, m_mink', u_k) - 0.5*SUMv;
S = H_x*P_mink*H_x' + 0.5*SUMs + R;
K = P_mink*H_x'/S;
m_k = m_mink + K*v;
P_k = P_mink - K*S*K';

if norm(P_k, "fro") > frobNormBound
    P_k = P_k./max(P_k, [], 'all').*rescaleFactor;
end

end

