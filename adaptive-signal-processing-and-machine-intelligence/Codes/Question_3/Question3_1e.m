fSample = 1e3;
fPower = 50;
nSamples = 1e3;
nPhases = 3;
amplitude = ones(nPhases, 1);
phaseShift = [0; -2/3 * pi; 2/3 * pi];
phaseDelay = [0; 0.2 * pi; 0.4 * pi];
ampDif = 0.4;
phaseInit = 0;
t = (0: nSamples - 1) / fSample;
orderFilter = 1;
step = 0.05;
leak = 0;

balancedABC = amplitude .* cos(2 * pi * fPower * t + phaseInit + phaseShift);
balancedZeroAlphaBeta = clarke(balancedABC);
balancedClarke = balancedZeroAlphaBeta(2, :) + 1i * balancedZeroAlphaBeta(3, :);
[circularityBalanced, ~] = circularity(balancedClarke);
[groupBalanced] = preprocessing(balancedClarke, orderFilter, 1);
[hBalancedClms, ~, errorBalancedClms] = clms(groupBalanced, balancedClarke, step, leak);
[hBalancedAclms, gBalancedAclms, ~, errorBalancedAclms] = aclms(groupBalanced, balancedClarke, step, leak);
fPowerBalancedClms = abs(fSample / (2 * pi) * atan(imag(hBalancedClms) ./ real(hBalancedClms)));
fPowerBalancedAclms = abs(fSample / (2 * pi) * atan(sqrt(imag(hBalancedAclms) .^2 - abs(gBalancedAclms) .^ 2) ./ real(hBalancedAclms)));

unbalancedPhaseABC = amplitude .* cos(2 * pi * fPower * t + phaseInit + phaseShift + phaseDelay);
unbalancedPhaseZeroAlphaBeta = clarke(unbalancedPhaseABC);
unbalancedPhaseClarke = unbalancedPhaseZeroAlphaBeta(2, :) + 1i * unbalancedPhaseZeroAlphaBeta(3, :);
[circularityUnbalancedPhase, ~] = circularity(unbalancedPhaseClarke);
[groupUnbalancedPhase] = preprocessing(unbalancedPhaseClarke, orderFilter, 1);
[hUnbalancedPhaseClms, ~, errorUnbalancedPhaseClms] = clms(groupUnbalancedPhase, unbalancedPhaseClarke, step, leak);
[hUnbalancedPhaseAclms, gUnbalancedPhaseAclms, ~, errorUnbalancedPhaseAclms] = aclms(groupUnbalancedPhase, unbalancedPhaseClarke, step, leak);
fPowerUnbalancedPhaseClms = abs(fSample / (2 * pi) * atan(imag(hUnbalancedPhaseClms) ./ real(hUnbalancedPhaseClms)));
fPowerUnbalancedPhaseAclms = abs(fSample / (2 * pi) * atan(sqrt(imag(hUnbalancedPhaseAclms) .^2 - abs(gUnbalancedPhaseAclms) .^ 2) ./ real(hUnbalancedPhaseAclms)));

unbalancedAmpABC = (amplitude + [-ampDif; 0; ampDif]) .* cos(2 * pi * fPower * t + phaseInit + phaseShift);
unbalancedAmpZeroAlphaBeta = clarke(unbalancedAmpABC);
unbalancedAmpClarke = unbalancedAmpZeroAlphaBeta(2, :) + 1i * unbalancedAmpZeroAlphaBeta(3, :);
[circularityUnbalancedAmp, ~] = circularity(unbalancedAmpClarke);
[groupUnbalancedAmp] = preprocessing(unbalancedAmpClarke, orderFilter, 1);
[hUnbalancedAmpClms, ~, errorUnbalancedAmpClms] = clms(groupUnbalancedAmp, unbalancedAmpClarke, step, leak);
[hUnbalancedAmpAclms, gUnbalancedAmpAclms, ~, errorUnbalancedAmpAclms] = aclms(groupUnbalancedAmp, unbalancedAmpClarke, step, leak);
fPowerUnbalancedAmpClms = abs(fSample / (2 * pi) * atan(imag(hUnbalancedAmpClms) ./ real(hUnbalancedAmpClms)));
fPowerUnbalancedAmpAclms = abs(fSample / (2 * pi) * atan(sqrt(imag(hUnbalancedAmpAclms) .^2 - abs(gUnbalancedAmpAclms) .^ 2) ./ real(hUnbalancedAmpAclms)));

figure;
subplot(3, 1, 1);
plot(fPowerBalancedClms, 'LineWidth', 2);
hold on;
plot(fPowerBalancedAclms, 'LineWidth', 2);
hold on;
plot([0 nSamples], [fPower fPower], 'k--', 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS', 'True');
title(sprintf('Frequency estimation for balanced system \\rho = 0'));
xlabel('Time (sample)');
ylabel('Frequency (Hz)');
ylim([0 100]);
subplot(3, 1, 2)
plot(fPowerUnbalancedPhaseClms, 'LineWidth', 2);
hold on;
plot(fPowerUnbalancedPhaseAclms, 'LineWidth', 2);
hold on;
plot([0 nSamples], [fPower fPower], 'k--', 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS', 'True');
title(sprintf('Frequency estimation for phase unbalanced system \\rho = %.2f', circularityUnbalancedPhase));
xlabel('Time (sample)');
ylabel('Frequency (Hz)');
ylim([0 100]);
subplot(3, 1, 3)
plot(fPowerUnbalancedAmpClms, 'LineWidth', 2);
hold on;
plot(fPowerUnbalancedAmpAclms, 'LineWidth', 2);
hold on;
plot([0 nSamples], [fPower fPower], 'k--', 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS', 'True');
title(sprintf('Frequency estimation for magnitude unbalanced system \\rho = %.2f', circularityUnbalancedAmp));
xlabel('Time (sample)');
ylabel('Frequency (Hz)');
ylim([0 100]);
figure;
subplot(3, 1, 1);
plot(pow2db(abs(errorBalancedClms) .^ 2), 'LineWidth', 2);
hold on;
plot(pow2db(abs(errorBalancedAclms) .^ 2), 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS');
title(sprintf('Learning curves for balanced system \\rho = 0'));
xlabel('Time (sample)');
ylabel('Error square (dB)');
subplot(3, 1, 2);
plot(pow2db(abs(errorUnbalancedPhaseClms) .^ 2), 'LineWidth', 2);
hold on;
plot(pow2db(abs(errorUnbalancedPhaseAclms) .^ 2), 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS');
title(sprintf('Learning curves for phase unbalanced system \\rho = %.2f', circularityUnbalancedPhase));
xlabel('Time (sample)');
ylabel('Error square (dB)');
subplot(3, 1, 3);
plot(pow2db(abs(errorUnbalancedAmpClms) .^ 2), 'LineWidth', 2);
hold on;
plot(pow2db(abs(errorUnbalancedAmpAclms) .^ 2), 'LineWidth', 2);
hold off;
grid on; grid minor;
legend('CLMS', 'ACLMS');
title(sprintf('Learning curves for magnitude unbalanced system \\rho = %.2f', circularityUnbalancedAmp));
xlabel('Time (sample)');
ylabel('Error square (dB)');