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

weightLms = cell(nSteps, nRps);
weightLmsAvg = cell(nSteps, 1);
for iStep = 1: nSteps
    for iRp = 1: nRps
        signal = arSignal(iRp, :);
        [group] = preprocessing(signal, orderAr, delay);
        [weightLms{iStep, iRp}, ~, ~] = lms(group, signal, step(iStep), leak);
    end
    weightLmsAvg{iStep} = mean(cat(3, weightLms{iStep, :}), 3);
end

legendStr = cell(2 * orderAr, 1);
figure;
for iStep = 1: nSteps
    subplot(nSteps, 1, iStep);
    for iOrder = 1: orderAr
        plot(weightLmsAvg{iStep}(iOrder, :), 'LineWidth', 1);
        hold on;
        legendStr{2 * iOrder - 1} = sprintf('Est. a_%d', iOrder);
        plot([0 nSamples], [coefAr(iOrder) coefAr(iOrder)], '--', 'LineWidth', 1);
        hold on;
        legendStr{2 * iOrder} = sprintf('a_%d', iOrder);
    end
    hold off;
    grid on; grid minor;
    legend(legendStr);
    title(sprintf('Steady state values of coefficients for \\mu = %.2f', step(iStep)));
    xlabel('Number of iterations');
    ylabel('Average weights');
    ylim([0 1]);
end