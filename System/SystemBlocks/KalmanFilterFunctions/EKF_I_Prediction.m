function [m_mink,P_mink] = EKF_I_Prediction(f, m_km, u_km, F_x, P_km, Q)
%EKF_I_PREDICTION Extended Kalman filter variant 1: first order with
%additive noise
%   Simo Särkkä (2013). Bayesian Filtering and Smoothing. Cambridge University Press. page 69
m_mink = f(m_km, u_km);
P_mink = F_x*P_km*F_x' + Q;
end

