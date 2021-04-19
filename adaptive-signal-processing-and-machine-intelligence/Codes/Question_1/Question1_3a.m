fSample = 1;
nSamples = 1024;
t = (0: nSamples - 1) / fSample;
fSine = [0.1 0.27];
noise = randn(1, nSamples);
noisySine = sin(2 * pi * fSine(1) * t) + sin(2 * pi * fSine(2) * t) + randn(1, nSamples);
filteredNoise = filter([1 1], 1, noise);
signal = {noise, noisySine, filteredNoise};

label = ["white Gaussian noise", "noisy sinusoid", "filtered white Gaussian noise"];
nSignals = length(signal);

acfUnbiased = cell(nSignals, 1);
acfBiased = cell(nSignals, 1);
psdUnbiased = cell(nSignals, 1);
psdBiased = cell(nSignals, 1);
for iSignal = 1: nSignals
    [acfUnbiased{iSignal}, lags] = xcorr(signal{iSignal}, 'unbiased');
    acfBiased{iSignal} = xcorr(signal{iSignal}, 'biased');
    psdUnbiased{iSignal} = real(fftshift(fft(ifftshift(acfUnbiased{iSignal}))));
    psdBiased{iSignal} = real(fftshift(fft(ifftshift(acfBiased{iSignal}))));
end
f = lags ./ (2 * nSamples) * fSample;


figure;
for iSignal = 1: nSignals
    subplot(nSignals, 1, iSignal);
    plot(lags, acfUnbiased{iSignal}, 'LineWidth', 2);
    hold on;
    plot(lags, acfBiased{iSignal}, 'LineWidth', 2);
    grid on; grid minor;
    legend('Unbiased', 'Biased');
    title(sprintf("Correlogram of %s", label(iSignal)));
    xlabel('Lag (sample)');
    ylabel('ACF');
end

figure;
for iSignal = 1: nSignals
    subplot(nSignals, 1, iSignal);
    plot(f, psdUnbiased{iSignal}, 'LineWidth', 2);
    hold on;
    plot(f, psdBiased{iSignal}, 'LineWidth', 2);
    grid on; grid minor;
    legend('Unbiased', 'Biased');
    title(sprintf("Spectral estimation by correlogram of %s", label(iSignal)));
    xlabel('Normalised frequency (\pi rad/sample)');
    ylabel('PSD');
end