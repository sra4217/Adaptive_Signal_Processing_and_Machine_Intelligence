ts = load('time-series.mat');
signal = (ts.y - mean(ts.y))';
step = 1e-5;
orderAr = 4;
delay = 1;
leak = 0;

[group] = preprocessing(signal, orderAr, delay);
[hLms, predictionLms, errorLms] = lms(group, signal, step, leak);
errorSquareLmsAvg = mean(abs(errorLms) .^ 2);
predGain = var(predictionLms) / var(errorLms);
predGainDb = pow2db(predGain);

figure;
plot(signal, 'k');
hold on;
plot(predictionLms, 'r');
hold off;
grid on; grid minor;
legend('Zero-mean', 'LMS');
title('Zero-mean signal and one-step ahead prediction by standard LMS');
xlabel('Time (sample)');
ylabel('Amplitude');
fprintf('MSE: %.4f dB\n', pow2db(errorSquareLmsAvg));
fprintf('Prediction gain %.4f dB\n', predGainDb);