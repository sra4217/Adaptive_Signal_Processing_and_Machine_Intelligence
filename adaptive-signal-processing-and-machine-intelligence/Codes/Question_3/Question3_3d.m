eeg = load('./data/EEG_Data_Assignment1.mat');
fSample = eeg.fs;
poz = detrend(eeg.POz - mean(eeg.POz))';
nSamples = 1200;
pos = 1e3;
poz = poz(pos: pos + nSamples - 1);
step = 1;
leak = [0, 1e-3];
nLeaks = length(leak);
dftMat = 1 / nSamples * exp(1i * (1: nSamples)' * pi / nSamples * (0: nSamples - 1));

psdDftClms = cell(nLeaks, 1);
for iLeak = 1: nLeaks
    [hFreqDftClms, ~, ~] = clms(dftMat, poz, step, leak(iLeak));
    psdDftClms{iLeak} = abs(hFreqDftClms) .^ 2;
    medianPsdDftClms = 1e3 * median(psdDftClms{iLeak}, 'all');
    psdDftClms{iLeak}(psdDftClms{iLeak} > medianPsdDftClms) = medianPsdDftClms;
end

figure;
for iLeak = 1: nLeaks
    subplot(nLeaks, 1, iLeak);
    mesh(psdDftClms{iLeak});
    view(2);
    cbar = colorbar;
    cbar.Label.String = 'PSD (dB)';
    grid on; grid minor;
    legend('DFT-CLMS');
    title([sprintf('Time-frequency diagram of EEG signal by DFT-CLMS \\mu = '), num2str(step), sprintf(' \\gamma = '), num2str(leak(iLeak))]);
    xlabel('Time (sample)');
    ylabel('Frequency (Hz)');
    ylim([0 70]);
end