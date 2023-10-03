function [IO, Llpmem] = footImpactDetector(Llpmem, stepcooldown)
    FDZeroCrossing = (Llpmem(2)-Llpmem(1))*(Llpmem(3)-Llpmem(2)) < 0; % first derivative zero crossing
%     ddLlp = ((Llpmem(3)-Llpmem(2))*120 - (Llpmem(2)-Llpmem(1))*120)*120;
%     IO = FDZeroCrossing && ddLlp < 0.005 && stepcooldown < 0;
    IO = FDZeroCrossing && stepcooldown < 0;
end