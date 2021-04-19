ts = load('time-series.mat');
signal = (ts.y - mean(ts.y))';
step = 1e-5;
orderAr = 4;
delay = 1;
leak = 0;
scale = 10: 10: 100;
nScales = length(scale);

predictionLms = cell(nScales, 1);
errorSquareLmsAvg = zeros(nScales, 1);
predGain = zeros(nScales, 1);
[group] = preprocessing(signal, orderAr, delay);
for iScale = 1: nScales
    [hLms, predictionLms{iScale}, errorLms] = lms_tanh(group, signal, step, leak, scale(iScale));
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
    legend('Zero-mean', 'Tanh-LMS');
    title(sprintf('One-step ahead prediction by tanh-LMS a = %d', scale(iScale)));
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
title('MSPE and prediction gain of tanh-LMS');