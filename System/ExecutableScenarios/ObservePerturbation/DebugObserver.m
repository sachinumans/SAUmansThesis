load UnperturbedRun.mat
xHatUnperturbed = xHat;
clearvars -except xHatUnperturbed
load PerturbedRun.mat
xHatPerturbed = xHat;
clearvars -except xHatUnperturbed xHatPerturbed
err = xHatUnperturbed - xHatPerturbed;

% close all;
figure()
for i = 1:3
    ax(i) = subplot(3,1,i);
    plot(err(i,:));
end
