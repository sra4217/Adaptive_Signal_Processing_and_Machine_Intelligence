ts = load('time-series.mat');
signal = ts.y';
nSamples = length(signal);
step = 1e-5;
orderAr = 4;
delay = 1;
leak = 0;
scale = 10: 10: 100;
nScales = length(scale);
nEpochs = 100;
nInits = 20;

hInit = cell(nScales, 1);
predictionLms = cell(nScales, 1);
errorSquareLmsAvg = zeros(nScales, 1);
predGain = zeros(nScales, 1);
batch = repmat(signal(1: nInits), 1, nEpochs);
[batchGroup] = preprocessing(batch, orderAr, delay);
augBatchGroup = [ones(1, size(batchGroup, 2)); batchGroup];
for iScale = 1: nScales
    [hInit{iScale}, ~, ~] = lms_tanh(augBatchGroup, batch, step, leak, scale(iScale));
    hInit{iScale} = hInit{iScale}(:, end);
end

[group] = preprocessing(signal, orderAr, delay);
augGroup = [ones(1, size(group, 2)); group];
for iScale = 1: nScales
    [~, predictionLms{iScale}, errorLms] = lms_tanh(augGroup, signal, step, leak, scale(iScale), hInit{iScale});
    errorSquareLmsAvg(iScale) = mean(abs(errorLms) .^ 2);
    predGain(iScale) = var(predictionLms{iScale}) / var(errorLms);
end
predGainDb = pow2db(predGain);

figure;
for iScale = 1: nScales
    subplot(nScales, 1, iScale);
    plot(signal, 'k');
    hold on;
    plot(predictionLms{iScale}, 'r');
    hold off;
    grid on; grid minor;
    legend('Non-zero-mean', 'Tanh-LMS');
    title(sprintf('One-step ahead prediction by biased tanh-LMS with pretrained weights a = %d', scale(iScale)));
    xlabel('Time (sample)');
    ylabel('Amplitude');
end
figure;
yyaxis left;
plot(scale, errorSquareLmsAvg, 'LineWidth', 2);
ylabel('MSPE (dB)');
yyaxis right;
plot(scale, predGainDb, 'LineWidth', 2);
ylabel('Prediction gain (dB)');
grid on; grid minor;
legend('MSPE', 'Prediction gain', 'location', 'northwest');
title('MSPE and prediction gain of pretrained biased tanh-LMS');
xlabel('Activation scale');