function [ResNorm] = compareModelPerStrideGA(p, p_bio, p_spring, w, k, xMeas, walkVel, gaitCycle, bound, LgrfPos, RgrfPos, LgrfVec, RgrfVec, LgrfMag, RgrfMag, LLML, LGTR, RLML, RGTR, dt, plotIO)
Wi = p_bio(1); l0 = p_bio(2);  m = p_bio(3); h = p_bio(4);
K_ss = p_spring(1);  b_ss = p_spring(2);  K_ds = p_spring(3);
Vl_ss = p(1); Vs_ss = p(2);
Vl_ds = p(3);
Vs_bl = p(4); Vs_fl = p(5);
J = p(6:8);
l_preload = p(9);

p = {0, 0, 0, 0, Wi, 0, h, l0+l_preload, [], [], [], [], [], []};
grfModel = cat(3, nan(3, 2, k(1)-1), nan(3, 2,length(k)));
legLenModel = [nan(2, k(1)-1), nan(2,length(k))];
k_switch = [];
linMult = [0];

ki = k(2);
xModel = [zeros(14, k(1)), xMeas(:,1), zeros(14,length(k)-2)];
while ki < k(end-1)
    k1 = ki;
    switch gaitCycle(1)
        case "lSS"
            k_end = ki+ find(RgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult (linMult(end)+1):(linMult(end)+k_end-ki)];
            lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            for ki = ki:k_end
                lFcur = lF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*LSSeom(0,xModel(:,ki)',lFcur,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                p{9} = K_ss; p{10} = b_ss;
                p{11} = Vs_ss; p{12} = Vl_ss;
                p{13} = [lF'; 0]; p{14} = 1;
                grfModel(:,1,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(1,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rSS"
            k_end = ki+ find(LgrfMag(ki:end)>bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult (linMult(end)+1):(linMult(end)+k_end-ki)];
            rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            for ki = ki:k_end
                rFcur = rF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*RSSeom(0,xModel(:,ki)',rFcur,Vl_ss,Vs_ss,h,Wi,l0+l_preload,m,K_ss,b_ss,J);

                p{9} = K_ss; p{10} = b_ss;
                p{11} = Vs_ss; p{12} = Vl_ss;
                p{13} = [rF'; 0]; p{14} = 2;
                grfModel(:,2,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(2,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "lDSr"
            k_end = ki+ find(LgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult 1:(k_end-ki)];
            lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");

            % reinitialise
            xModel(:, ki) = xMeas(:,ki-k(1)+1);
            for ki = ki:k_end
                lFcur = lF-walkVel(1:2)*dt*(ki-k1);
                rFcur = rF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*lDSr_split_eom(0,xModel(:,ki)',lFcur,rFcur,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,0,J);

                p{9} = K_ds; p{10} = 0;
                p{11} = [Vs_bl, Vs_fl]; p{12} = Vl_ds;
                p{13} = [[lF'; 0],[rF'; 0]]; p{14} = -1;
                grfModel(:,:,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);
        case "rDSl"
            k_end = ki+ find(RgrfMag(ki:end) < bound, 1);
            k_switch = [k_switch k_end];
            linMult = [linMult 1:(k_end-ki)];
            lF = mean(LgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");
            rF = mean(RgrfPos(k1:k_end,1:2) + (walkVel(1:2)'*(0:(k_end-k1)))'*dt, 1, "omitnan");

            % reinitialise
            xModel(:, ki) = xMeas(:,ki-k(1)+1);
            for ki = ki:k_end
                lFcur = lF-walkVel(1:2)*dt*(ki-k1);
                rFcur = rF-walkVel(1:2)*dt*(ki-k1);
                xModel(:,ki+1) = xModel(:,ki) + dt*rDSl_split_eom(0,xModel(:,ki)',lFcur,rFcur,Vl_ds,Vs_bl,Vs_fl,h,Wi,l0+l_preload,m,K_ds,0,J);

                p{9} = K_ds; p{10} = 0;
                p{11} = [Vs_bl, Vs_fl]; p{12} = Vl_ds;
                p{13} = [[lF'; 0],[rF'; 0]]; p{14} = -2;
                grfModel(:,:,ki) = state2grf_DSsplit(xModel(:,ki), p);

                legLenModel(:,ki) = state2legLength_DSsplit(xModel(:,ki), p);
            end
            gaitCycle = circshift(gaitCycle, -1);

    end
end

xModel = xModel(:,k(2:end-1));
xModelRes = w*(xModel - xMeas)*diag(linMult(2:(length(xModel)+1))'.^2);
xModelRes(isnan(xModelRes)) = 1e7;
xModelResNorm = norm(xModelRes, "fro")^2;


ResNorm = xModelResNorm;

if plotIO
    plotModelMeasComparison(k, dt, grfModel, LgrfVec, RgrfVec, legLenModel,...
        LLML, LGTR, RLML, RGTR, xMeas, xModel, k_switch, gaitCycle);
end

end