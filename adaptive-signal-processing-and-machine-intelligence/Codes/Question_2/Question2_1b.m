nSamples = 1e3;
nRps = 1e2;
coefAr = [0.1 0.8];
orderAr = length(coefAr);
variance = 0.25;
delay = 1;
step = [0.05; 0.01];
nSteps = length(step);
leak = 0;

arModel = arima('AR', coefAr, 'Variance', variance, 'Constant', 0);
arSignal = simulate(arModel, nSamples, 'NumPaths', nRps);
arSignal = arSignal';

error = cell(nSteps, nRps);
errorSquareAvg = cell(nSteps, 1);
for iStep = 1: nSteps
    for iRp = 1: nRps
        signal = arSignal(iRp, :);
        [group] = preprocessing(signal, orderAr, delay);
        [~, ~, error{iStep, iRp}] = lms(group, signal, step(iStep), leak);
    end
    errorSquareAvg{iStep} = mean(cat(3, error{iStep, :}) .^ 2, 3);
end

legendStr = cell(nSteps, 1);
figure;
subplot(2, 1, 1);
for iStep = 1: nSteps
    plot(pow2db(error{iStep, end}.^2), 'LineWidth', 1);
    legendStr{iStep} = sprintf('\\mu = %.2f', step(iStep));
    hold on;
end
hold off;
grid on; grid minor;
legend(legendStr, 'location', 'southeast');
title('Error instance by adaptive LMS with second-order AR model');
xlabel('Time');
ylabel('Squared Error');
subplot(2, 1, 2);
for iStep = 1: nSteps
    plot(pow2db(errorSquareAvg{iStep}), 'LineWidth', 1);
    legendStr{iStep} = sprintf('\\mu = %.2f', step(iStep));
    hold on;
end
hold off;
grid on; grid minor;
legend(legendStr, 'location', 'northeast');
title('Mean error by adaptive LMS with second-order AR model');
xlabel('Time');
ylabel('Squared Error');