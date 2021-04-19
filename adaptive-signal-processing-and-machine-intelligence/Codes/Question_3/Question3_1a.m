nSamples = 1e3;
nRps = 1e2;
coefMa = [1.5+1i, 2.5-0.5i];
orderFilter = 2;
step = 0.1;
leak = 0;

whiteNoise = 1 / sqrt(2) * (randn(nRps, nSamples) + 1i * randn(nRps, nSamples));
[circularityCoefWn, ~] = circularity(whiteNoise);
wlmaSignal = whiteNoise + coefMa(1) * [zeros(nRps, 1), whiteNoise(:, 1: end - 1)] + coefMa(2) * [zeros(nRps, 1), conj(whiteNoise(:, 1: end - 1))];
[circularityCoefWlma, ~] = circularity(wlmaSignal);

errorClms = cell(nRps, 1);
errorAclms = cell(nRps, 1);
for iRp = 1: nRps
    delayedSignal = [0, (wlmaSignal(iRp, 1: end - 1))];
    [group] = preprocessing(whiteNoise(iRp, :), orderFilter, 1);
    [~, ~, errorClms{iRp}] = clms(group, delayedSignal, step, leak);
    [~, ~, ~, errorAclms{iRp}] = aclms(group, delayedSignal, step, leak);
end
errorSquareClmsAvg = mean(abs(cat(3, errorClms{:})) .^ 2, 3);
errorSquareAclmsAvg = mean(abs(cat(3, errorAclms{:})) .^ 2, 3);

figure;
subplot(2, 2, 1);
scatter(real(whiteNoise(:)), imag(whiteNoise(:)), 1, 'k');
legend('White noise');
title(sprintf('White noise \\rho = %.2f', circularityCoefWn));
xlabel('Real part');
ylabel('Imaginary part');
subplot(2, 2, 2);
scatter(real(wlmaSignal(:)), imag(wlmaSignal(:)), 1, 'r');
legend('WLMA (1)');
title(sprintf('WLMA (1) \\rho = %.2f', circularityCoefWlma));
xlabel('Real part');
ylabel('Imaginary part');
subplot(2, 1, 2);
plot(pow2db(errorSquareClmsAvg), 'LineWidth', 2);
hold on;
plot(pow2db(errorSquareAclmsAvg), 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS');
title('Error curves of CLMS and ACLMS');
xlabel('Time (sample)');
ylabel('Squared error (dB)');
