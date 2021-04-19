fSample = 1;
nSamples = 1024;
t = (0: nSamples - 1) / fSample;
f = (-nSamples / 2: nSamples / 2 - 1) * (fSample / nSamples);
fExp1 = 0.3;
expWave{1} = exp(1i * 2 * pi * fExp1 * t);
fExp2 = [0.3, 0.32];
expWave{2} = exp(1i * 2 * pi * fExp2(1) * t) + exp(1i * 2 * pi * fExp2(2) * t);
fExp3 = [0.28, 0.3, 0.32];
expWave{3} = exp(1i * 2 * pi * fExp3(1) * t) + exp(1i * 2 * pi * fExp3(2) * t) + exp(1i * 2 * pi * fExp3(3) * t);
nExps = 3;
point = 20: 10: 50;
nPoints = length(point);
pNoise = 0.2;


psd = cell(nExps, nPoints);
for iExp = 1: nExps
    for iPoint = 1: nPoints
        noise = sqrt(pNoise / 2) * (randn(1, point(iPoint)) + 1i * randn(1, point(iPoint)));
        noisyExp = expWave{iExp}(1: point(iPoint)) + noise;
        psd{iExp, iPoint} = abs(fftshift(fft(noisyExp, nSamples))) .^ 2 / point(iPoint);
    end
end


legendStr = cell(nPoints, 1);
figure;
for iExp = 1: nExps
    subplot(nExps, 1, iExp);
    for iPoint = 1: nPoints
        plot(f, psd{iExp, iPoint}, 'LineWidth', 2);
        legendStr{iPoint} = ['N = ', num2str(point(iPoint))];
        hold on;
    end
    grid on; grid minor;
    title(sprintf('Periodogram of signal with %d exponential(s)', iExp));
    legend(legendStr);
    xlabel('Normalised frequency');
    ylabel('PSD');
    xlim([0.25, 0.4]);
end
