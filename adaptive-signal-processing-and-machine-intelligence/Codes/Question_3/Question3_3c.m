fSample = 1500;
nSamples = 1500;
variance = 0.05;
nSegs = 3;
lengthSeg = 500;
freqFunc = @(n) ((1 <= n) & (n <= 500)) .* 100 + ((501 <= n) & (n <= 1000)) .* (100 + (n - 500) / 2) + ((1001 <= n) & (n <= 1500)) .* (100 + ((n - 1000) / 25) .^ 2);
freqSeq = freqFunc(1: nSamples);
phaseSeq = cumsum(freqSeq);
orderFilter = 1;
stepAr = 0.1; stepDft = 1;
leak = [0, 1e-2, 1e-1, 5e-1];
nLeaks = length(leak);
fmSignal = exp(1i * 2 * pi / fSample * phaseSeq) + sqrt(variance / 2) * (randn(1, nSamples) + 1i * randn(1, nSamples));
dftMat = 1 / nSamples * exp(1i * (1: nSamples)' * 2 * pi / nSamples * (0: nSamples - 1));

psdArClms = cell(nLeaks, 1);
[group] = preprocessing(fmSignal, orderFilter, 1);
for iLeak = 1: nLeaks
    [hArClms, ~, ~] = clms(group, fmSignal, stepAr, leak(iLeak));
    for iSample = 1: nSamples
        [hFreqArClms, fArClms] = freqz(1, [1; -conj(hArClms(iSample))], nSamples, fSample);
        psdArClms{iLeak}(:, iSample) = abs(hFreqArClms) .^ 2;
    end
    medianPsdArClms = 50 * median(psdArClms{iLeak}, 'all');
    psdArClms{iLeak}(psdArClms{iLeak} > medianPsdArClms) = medianPsdArClms;
end

psdDftClms = cell(nLeaks, 1);
for iLeak = 1: nLeaks
    [hFreqDftClms, ~, ~] = clms(dftMat, fmSignal, stepDft, leak(iLeak));
    psdDftClms{iLeak} = abs(hFreqDftClms) .^ 2;
    medianPsdDftClms = 50 * median(psdDftClms{iLeak}, 'all');
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
    title([sprintf('Time-frequency diagram of FM signal by DFT-CLMS \\mu = '), num2str(stepDft), sprintf(' \\gamma = '), num2str(leak(iLeak))]);
    xlabel('Time (sample)');
    ylabel('Frequency (Hz)');
    ylim([0 1000]);
end