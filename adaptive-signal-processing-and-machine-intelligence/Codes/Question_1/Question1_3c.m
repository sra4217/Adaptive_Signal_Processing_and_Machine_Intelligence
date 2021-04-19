fSample = 1;
nSamples = 1024;
t = (0: nSamples - 1) / fSample;
aSine = [0.7 0.5];
fSine = [0.1 0.27];
sineWave = aSine(1) * sin(2 * pi * fSine(1) * t) + aSine(2) * sin(2 * pi * fSine(2) * t);
nRps = 1e2;

acf = cell(nRps, 1);
psd = cell(nRps, 1);
for iRp = 1: nRps
    noisySine = sineWave + randn(size(sineWave));
    [acf{iRp}, lags] = xcorr(noisySine, 'biased');
    psd{iRp} = real(fftshift(fft(ifftshift(acf{iRp}))));
end
f = lags ./ (2 * nSamples) * fSample;
psdMean = mean(cell2mat(psd));
psdStd = std(cell2mat(psd));

figure;
subplot(2, 1, 1);
for iRp = 1: nRps
    irPlot = plot(f, pow2db(psd{iRp}), 'k', 'LineWidth', 2);
    hold on;
end
meanPlot = plot(f, pow2db(psdMean), 'LineWidth', 2);
hold off;
grid on; grid minor;
legend([irPlot, meanPlot], {'Individual', 'Mean'}, 'location', 'southeast');
title('Individual and mean PSD using biased estimator');
xlabel('Normalised frequency');
ylabel('PSD');
subplot(2, 1, 2);
varPlot = plot(f, pow2db(psdStd), 'r', 'LineWidth', 2);
grid on; grid minor;
legend('Standard deviation');
title('Standard deviation of PSD estimate of noise-corrupted signals');
xlabel('Normalised frequency');
ylabel('PSD');
