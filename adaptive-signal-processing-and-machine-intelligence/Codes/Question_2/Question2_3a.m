fSample = 1;
nSamples = 1e3;
t = (0: nSamples - 1) / fSample;
aSine = 1;
fSine = 5e-3;
signal = aSine * sin(2 * pi * fSine * t);
nRps = 1e2;
coefMa = [0 0.5];
variance = 1;
step = 0.01;
nDelays = 7;
orderFilter = 5;
leak = 0;
nTransients = 50;

maModel = arima('MA', coefMa, 'Variance', variance, 'Constant', 0);
[maSignal, innovation] = simulate(maModel, nSamples, 'NumPaths', nRps);
colouredNoise = maSignal';
whiteNoise = innovation';

noisySignal = cell(nDelays, nRps);
signalAle = cell(nDelays, nRps);
errorSquare = cell(nDelays, nRps);
mspe = zeros(nDelays, 1);
for iDelay = 1: nDelays
    for iRp = 1: nRps
        noisySignal{iDelay, iRp} = signal + colouredNoise(iRp, :);
        [group] = preprocessing(noisySignal{iDelay, iRp}, orderFilter, iDelay);
        [~, signalAle{iDelay, iRp}, ~] = lms(group, noisySignal{iDelay, iRp}, step, leak);
        errorSquare{iDelay, iRp} = (signal(nTransients + 1: end) - signalAle{iDelay, iRp}(nTransients + 1: end)) .^ 2;
    end
    mspe(iDelay) = mean(cell2mat(errorSquare(iDelay, :)));
end

figure;
for iDelay = 1: nDelays
    subplot(nDelays, 1, iDelay);
    for iRp = 1: nRps
        noisyPlot = plot(t, noisySignal{iDelay, iRp}, 'k', 'LineWidth', 2);
        hold on;
        alePlot = plot(t, signalAle{iDelay, iRp}, 'b', 'LineWidth', 2);
        hold on;
    end
    cleanPlot = plot(t, signal, 'r', 'LineWidth', 2);
    hold off;
    grid on; grid minor;
    legend([noisyPlot, alePlot, cleanPlot], {'Noisy', 'ALE', 'Clean'}, 'location', 'bestoutside');
    title(sprintf('Noisy, clean and ALE signals by linear predictor M = %d \\Delta = %d', orderFilter, iDelay));
    xlabel('Time (sample)');
    ylabel('Amplitude');
end
figure;
plot(pow2db(mspe), 'm', 'LineWidth', 2);
grid on; grid minor;
legend('MSPE');
title(sprintf('MSPE by linear predictor M = %d', orderFilter));
xlabel('Delay (sample)');
ylabel('MSPE (dB)');