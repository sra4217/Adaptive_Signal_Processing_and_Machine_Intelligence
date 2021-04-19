eeg = load('./data/EEG_Data_Assignment2.mat');
fSample = eeg.fs;
nSamples = length(eeg.POz);
nFft = 2 ^ 13;
poz = (eeg.POz - mean(eeg.POz))';
variance = 0.01;
t = (0: nSamples - 1) / fSample;
tWindow = 5;
nWindows = 2 ^ 12;
nOverlaps = round(0.9 * nWindows);
step = [1e-2, 1e-3, 1e-4];
nSteps = length(step);
orderFilter = 5: 5: 15;
nOrders = length(orderFilter);
leak = 0;
nTransients = 50;

ampSine = 1;
freqSine = 50;
sine = ampSine * sin(2 * pi * freqSine * t);
noise = variance * randn(1, nSamples) + sine;

signalAnc = cell(nOrders, nSteps);
mspeAnc = zeros(nOrders, nSteps);
errorSquareAnc = cell(nOrders, nSteps);
delayedPoz = [0, poz(1: end - 1)];
for iOrder = 1: nOrders
    [group] = preprocessing(noise, orderFilter(iOrder), 1);
    for iStep = 1: nSteps
        [~, noiseAnc, ~] = lms(group, delayedPoz, step(iStep), leak);
        signalAnc{iOrder, iStep} = delayedPoz - noiseAnc;
        errorSquareAnc{iOrder, iStep} = (poz(nTransients + 1: end) - signalAnc{iOrder, iStep}(nTransients + 1: end)) .^ 2;
        mspeAnc(iOrder, iStep) = mean(errorSquareAnc{iOrder, iStep});
    end
end
orderOptimal = 10;
stepOptimal = 1e-3;
orderIndex = find(orderFilter == orderOptimal);
stepIndex = find(step == stepOptimal);
[psdPoz, fAnalog] = periodogram(poz, rectwin(nSamples), nFft, fSample);
psdAnc = periodogram(signalAnc{orderIndex, stepIndex}, rectwin(nSamples), nFft, fSample);


figure;
spectrogram(poz, nWindows, nOverlaps, nFft, fSample, 'yaxis');
ylim([0 60]);
title('Spectrogram of preprocessed POz');
for iOrder = 1: nOrders
    figure;
    for iStep = 1: nSteps
        subplot(nSteps, 1, iStep);
        spectrogram(signalAnc{iOrder, iStep}, nWindows, nOverlaps, nFft, fSample, 'yaxis');
        ylim([0 60]);
        title(['Spectrogram of ANC signal by linear predictor M = ', num2str(orderFilter(iOrder)), sprintf(' and \\mu = '), num2str(step(iStep))]);
    end
end
figure;
legendStr = cell(nSteps, 1);
for iStep = 1: nSteps
    plot(orderFilter, pow2db(mspeAnc(:, iStep)), 'LineWidth', 2);
    legendStr{iStep} = [sprintf('\\mu = '), num2str(step(iStep))];
    hold on;
end
hold off;
grid on; grid minor;
legend(legendStr, 'location', 'southeast');
title('MSPE against filter order and step size');
xlabel('Order');
ylabel('MSPE (dB)');
figure;
subplot(2, 1, 1);
plot(fAnalog, pow2db(psdPoz), 'LineWidth', 2);
hold on;
plot(fAnalog, pow2db(psdAnc), 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('Original', 'ANC');
title('Periodograms of original and optimal ANC POz');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
xlim([0 60]);
ylim([-160 -100]);
subplot(2, 1, 2);
plot(fAnalog, pow2db(abs(psdPoz - psdAnc)), 'm', 'LineWidth', 2);
grid on; grid minor;
legend('Absolute error');
title('Periodogram error by optimal ANC');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
xlim([0 60]);
ylim([-160 -100]);