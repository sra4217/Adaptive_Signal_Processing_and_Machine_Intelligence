fSample = 1;
nSamples = 1e3;
t = (0: nSamples - 1) / fSample;
aSine = 1;
fSine = 5e-3;
signal = aSine * sin(2 * pi * fSine * t);
nRps = 1e2;
coefMa = [0 0.5];
variance = 1;
step = 0.005;
delay = 3;
nDelays = length(delay);
orderFilter = 6;
leak = 0;
nTransients = 50;


maModel = arima('MA', coefMa, 'Variance', variance, 'Constant', 0);
[maSignal, innovation] = simulate(maModel, nSamples, 'NumPaths', nRps);
colouredNoise = maSignal';
whiteNoise = innovation';
secondaryNoise = 0.9 * colouredNoise + 0.05;

noisySignal = cell(nDelays, nRps);
signalAle = cell(nDelays, nRps);
signalAnc = cell(nDelays, nRps);
errorSquareAle = cell(nDelays, nRps);
errorSquareAnc = cell(nDelays, nRps);
mspeAle = zeros(nDelays, 1);
mspeAnc = zeros(nDelays, 1);
signalAleAvg = cell(nDelays, 1);
signalAncAvg = cell(nDelays, 1);
for iDelay = 1: nDelays
    for iRp = 1: nRps
        noisySignal{iDelay, iRp} = signal + colouredNoise(iRp, :);
        [group] = preprocessing(noisySignal{iDelay, iRp}, orderFilter, delay(iDelay));
        [~, signalAle{iDelay, iRp}, ~] = lms(group, noisySignal{iDelay, iRp}, step, leak);
        errorSquareAle{iDelay, iRp} = (signal(nTransients + 1: end) - signalAle{iDelay, iRp}(nTransients + 1: end)) .^ 2;
        delayedSignal = [0, noisySignal{iDelay, iRp}(1: end - 1)];
        [group] = preprocessing(secondaryNoise(iRp, :), orderFilter, 1);
        [~, noiseAnc, ~] = lms(group, delayedSignal, step, leak);
        signalAnc{iDelay, iRp} = delayedSignal - noiseAnc;
        errorSquareAnc{iDelay, iRp} = (signal(nTransients + 1: end) - signalAnc{iDelay, iRp}(nTransients + 1: end)) .^ 2;
    end
    mspeAle(iDelay) = mean(cell2mat(errorSquareAle(iDelay, :)));
    mspeAnc(iDelay) = mean(cell2mat(errorSquareAnc(iDelay, :)));
    signalAleAvg{iDelay} = mean(cat(3, signalAle{iDelay, :}), 3);
    signalAncAvg{iDelay} = mean(cat(3, signalAnc{iDelay, :}), 3);
end


figure;
subplot(3, 1, 1);
for iDelay = 1: nDelays
    for iRp = 1: nRps
        noisyPlot = plot(t, noisySignal{iDelay, iRp}, 'k', 'LineWidth', 2);
        hold on;
        alePlot = plot(t, signalAle{iDelay, iRp}, 'b', 'LineWidth', 2);
        hold on;
    end
    cleanPlot = plot(t, signal, 'r', 'LineWidth', 2);
    hold off;
    grid on; grid minor;
    legend([noisyPlot, alePlot, cleanPlot], {'Noisy', 'ALE', 'Clean'});
    title(sprintf('ALE signal by linear predictor M = %d and \\Delta = %d, MSPE = %.2f dB', orderFilter, delay(iDelay), pow2db(mspeAle)));
    xlabel('Time (sample)');
    ylabel('Amplitude');
end
subplot(3, 1, 2);
for iDelay = 1: nDelays
    for iRp = 1: nRps
        noisyPlot = plot(t, noisySignal{iDelay, iRp}, 'k', 'LineWidth', 2);
        hold on;
        ancPlot = plot(t, signalAnc{iDelay, iRp}, 'c', 'LineWidth', 2);
        hold on;
    end
    cleanPlot = plot(t, signal, 'r', 'LineWidth', 2);
    hold off;
    grid on; grid minor;
    legend([noisyPlot, ancPlot, cleanPlot], {'Noisy', 'ANC', 'Clean'});
    title(sprintf('ANC signal by linear predictor M = %d and \\Delta = %d, MSPE = %.2f dB', orderFilter, delay(iDelay), pow2db(mspeAnc)));
    xlabel('Time (sample)');
    ylabel('Amplitude');
end
subplot(3, 1, 3);
for iDelay = 1: nDelays
    hold on;
    alePlot = plot(t, signalAleAvg{iDelay}, 'LineWidth', 2);
    hold on;
    ancPlot = plot(t, signalAncAvg{iDelay}, 'LineWidth', 2);
    hold on;
    cleanPlot = plot(t, signal, 'k', 'LineWidth', 2);
    hold off;
    grid on; grid minor;
    legend([alePlot, ancPlot, cleanPlot], {'ALE', 'ANC', 'Clean'});
    title(sprintf('ALE, ANC and clean signals by linear predictor M = %d and \\Delta = %d', orderFilter, delay(iDelay)));
    xlabel('Time (sample)');
    ylabel('Amplitude');
end
