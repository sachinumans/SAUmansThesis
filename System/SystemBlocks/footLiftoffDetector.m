function [IO] = footLiftoffDetector(Llp_mem, liftcooldown)
    ZeroCrossing = (Llp_mem(3)*Llp_mem(2)) < 0;
    dLlp = (Llp_mem(3)-Llp_mem(2))*120;
    IO = ZeroCrossing && dLlp > -0.005 && liftcooldown < 0;
end