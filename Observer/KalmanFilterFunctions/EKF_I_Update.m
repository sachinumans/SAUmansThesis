function [m_k, P_k] = EKF_I_Update(y, h, m_mink, u_k, H_x, P_mink, R, frobNormBound, rescaleFactor)
%EKF_I_UPDATE Extended Kalman filter variant 1: first order with
%additive noise
%   Simo Särkkä (2013). Bayesian Filtering and Smoothing. Cambridge University Press. page 70
v = y - h(0, m_mink', u_k);
S = H_x*P_mink*H_x' + R;
K = P_mink*H_x'/S;
m_k = m_mink + K*v;
P_k = P_mink - K*S*K';

if norm(P_k, "fro") > frobNormBound % reset covariance
    P_k = P_k./max(P_k, [], 'all').*rescaleFactor;
end

end

