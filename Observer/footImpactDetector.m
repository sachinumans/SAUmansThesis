function [IO, Llpmem] = footImpactDetector(fi, Llpmem, stepcooldown)
    Llpmem = circshift(Llpmem, -1); Llpmem(3) = fi(end);

    FDZeroCrossing = (Llpmem(2)-Llpmem(1))*(Llpmem(3)-Llpmem(2)) < 0;
    ddLlp = ((Llpmem(3)-Llpmem(2))*120 - (Llpmem(2)-Llpmem(1))*120)*120;
    IO = FDZeroCrossing && ddLlp > -0.005 && stepcooldown < 0;
end