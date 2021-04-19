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
delay = 3: 25;
nDelays = length(delay);
orderFilter = 5: 5: 20;
nOrders = length(orderFilter);
leak = 0;
nTransients = 50;

maModel = arima('MA', coefMa, 'Variance', variance, 'Constant', 0);
[maSignal, innovation] = simulate(maModel, nSamples, 'NumPaths', nRps);
colouredNoise = maSignal';
whiteNoise = innovation';

errorSquare = cell(nOrders, nDelays, nRps);
mspe = zeros(nOrders, nDelays);
for iOrder = 1: nOrders
    for iDelay = 1: nDelays
        for iRp = 1: nRps
            noisySignal = signal + colouredNoise(iRp, :);
            [group] = preprocessing(noisySignal, orderFilter(iOrder), delay(iDelay));
            [~, signalAle, ~] = lms(group, noisySignal, step, leak);
            errorSquare{iOrder, iDelay, iRp} = (signal(nTransients + 1: end) - signalAle(nTransients + 1: end)) .^ 2;
        end
        mspe(iOrder, iDelay) = mean(cell2mat(errorSquare(iOrder, iDelay, :)), 'all');
    end
end

figure;
legendStr = cell(nOrders, 1);
subplot(2, 1, 1);
for iOrder = 1: nOrders
    plot(delay, pow2db(mspe(iOrder, :)), 'LineWidth', 2);
    legendStr{iOrder} = sprintf('M = %d', orderFilter(iOrder));
    hold on;
end
grid on; grid minor;
legend(legendStr, 'location', 'southeast');
title('MSPE against delay and filter order');
xlabel('Delay (sample)');
ylabel('MSPE (dB)');
xlim([min(delay), max(delay)]);
subplot(2, 1, 2);
nDelayPlots = 1;
legendStr = cell(nDelayPlots, 1);
for iDelayPlot = 1: nDelayPlots
    plot(orderFilter, pow2db(mspe(:, iDelayPlot)), 'LineWidth', 2);
    legendStr{iDelayPlot} = sprintf('\\Delta = %d', delay(iDelayPlot));
    hold on;
end
grid on; grid minor;
legend(legendStr, 'location', 'southeast');
title('MSPE against filter order');
xlabel('Filter order');
ylabel('MSPE (dB)');
