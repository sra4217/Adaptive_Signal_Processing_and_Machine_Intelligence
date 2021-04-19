fSample = 1500;
nSamples = 1500;
variance = 0.05;
nSegs = 3;
lengthSeg = 500;
freqFunc = @(n) ((1 <= n) & (n <= 500)) .* 100 + ((501 <= n) & (n <= 1000)) .* (100 + (n - 500) / 2) + ((1001 <= n) & (n <= 1500)) .* (100 + ((n - 1000) / 25) .^ 2);
freqSeq = freqFunc(1: nSamples);
phaseSeq = cumsum(freqSeq);
orderFilter = 1;
step = [1e0, 1e-1, 1e-2, 1e-3];
nSteps = length(step);
leak = 0;
fmSignal = exp(1i * 2 * pi / fSample * phaseSeq) + sqrt(variance / 2) * (randn(1, nSamples) + 1i * randn(1, nSamples));

psdArClms = cell(nSteps, 1);
[group] = preprocessing(fmSignal, orderFilter, 1);
for iStep = 1: nSteps
    [hArClms, ~, ~] = clms(group, fmSignal, step(iStep), leak);
    for iSample = 1: nSamples
        [hFreqArClms, fArClms] = freqz(1, [1; -conj(hArClms(iSample))], nSamples, fSample);
        psdArClms{iStep}(:, iSample) = abs(hFreqArClms) .^ 2;
    end
    medianPsdArClms = 50 * median(psdArClms{iStep}, 'all');
    psdArClms{iStep}(psdArClms{iStep} > medianPsdArClms) = medianPsdArClms;
end

figure;
for iStep = 1: nSteps
    subplot(nSteps, 1, iStep);
    mesh(psdArClms{iStep});
    view(2);
    cbar = colorbar;
    cbar.Label.String = 'PSD (dB)';
    grid on; grid minor;
    legend(sprintf('CLMS-AR (%d)', orderFilter));
    title([sprintf('Time-frequency diagram of FM signal by CLMS \\mu = '), num2str(step(iStep))]);
    xlabel('Time (sample)');
    ylabel('Frequency (Hz)');
end
