ecg = load('./data/ECG_Data_Assignment1.mat');
fSample = ecg.fsRRI;
rri = {detrend(ecg.xRRI1 - mean(ecg.xRRI1)) detrend(ecg.xRRI2 - mean(ecg.xRRI2)) detrend(ecg.xRRI3 - mean(ecg.xRRI3))};
nRris = length(rri);
label = ["normal", "fast", "slow"];
tWindow = [50 150];
nFft = 2048;
nOverlap = 0;

psdStd = cell(nRris, 1);
for iRri = 1: nRris
    nSamples = length(rri{iRri});
    [psdStd{iRri}, f] = periodogram(rri{iRri}, hamming(nSamples), nFft, fSample);
end

psdAvg = cell(nRris, length(tWindow));
for iRri = 1: nRris
    for iWindow = 1: length(tWindow)
        nWindows = tWindow(iWindow) * fSample;
        psdAvg{iRri, iWindow} = pwelch(rri{iRri}, hamming(nWindows), nOverlap, nFft, fSample);
    end
end

legendStr = cell(1, length(tWindow) + 1);
figure;
for iRri = 1: nRris
    subplot(nRris, 1, iRri);
    plot(f, pow2db(psdStd{iRri}), 'LineWidth', 2);
    hold on;
    legendStr{1} = 'Standard';
    for iWindow = 1: length(tWindow)
        plot(f, pow2db(psdAvg{iRri, iWindow}), 'LineWidth', 1);
        hold on;
        legendStr{iWindow + 1} = sprintf('\\Deltat = %d sec', tWindow(iWindow));
    end
    grid on; grid minor;
    legend(legendStr);
    title(sprintf('Periodogram by standard and Bartlett methods for %s RRI', label(iRri)));
    xlabel('Frequency');
    ylabel('PSD');
    ylim([-80 0]);
end
