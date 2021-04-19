load sunspot.dat
fSample = 1;
t = sunspot(:, 1);
xOriginal = sunspot(:, 2);
nSamples = length(t);
nOverlap = 0;

% removing mean and detrend
xMeanDetrend = detrend(xOriginal - mean(xOriginal));
% applying log then subtract mean
xLogMean = log(xOriginal + eps) - mean(log(xOriginal + eps));

[psdRaw, ~] = pwelch(xOriginal, hamming(nSamples), nOverlap, nSamples, fSample);
[psdMeanDetrend, ~] = pwelch(xMeanDetrend, hamming(nSamples), nOverlap, nSamples, fSample);
[psdLogMean, f] = pwelch(xLogMean, hamming(nSamples), nOverlap, nSamples, fSample);

figure;
subplot(2, 1, 1);
plot(t, xOriginal, 'LineWidth', 2);
hold on;
plot(t, xMeanDetrend, 'LineWidth', 2);
hold on;
plot(t, xLogMean, 'LineWidth', 2);
grid on; grid minor;
legend('Original', 'Mean-detrend', 'Log-mean');
title('Sunspot time series');
xlabel('Year');
ylabel('Number of sunspots');

subplot(2, 1, 2);
plot(f, pow2db(psdRaw), 'LineWidth', 2);
hold on;
plot(f, pow2db(psdMeanDetrend), '--', 'LineWidth', 3);
hold on;
plot(f, pow2db(psdLogMean), 'LineWidth', 2);
grid on; grid minor;
legend('Original', 'Mean-detrend', 'Log-mean');
title('Sunspots periodogram with Hamming window');
xlabel('Normalised frequency');
ylabel('PSD');
ylim([-30 70]);