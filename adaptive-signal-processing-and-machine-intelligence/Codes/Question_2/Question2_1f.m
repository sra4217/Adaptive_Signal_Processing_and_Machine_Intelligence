nSamples = 1e3;
nRps = 1e2;
coefAr = [0.1 0.8];
orderAr = length(coefAr);
variance = 0.25;
delay = 1;
step = [0.05; 0.01];
nSteps = length(step);
leak = 0.2: 0.3: 0.8;
nLeaks = length(leak);

arModel = arima('AR', coefAr, 'Variance', variance, 'Constant', 0);
arSignal = simulate(arModel, nSamples, 'NumPaths', nRps);
arSignal = arSignal';

weightLeakyLms = cell(nLeaks, nSteps, nRps);
weightLeakyLmsAvg = cell(nLeaks, nSteps);
for iLeak = 1: nLeaks
    for iStep = 1: nSteps
        for iRp = 1: nRps
            signal = arSignal(iRp, :);
            [group] = preprocessing(signal, orderAr, delay);
            [weightLeakyLms{iLeak, iStep, iRp}, ~, ~] = lms(group, signal, step(iStep), leak(iLeak));
        end
        weightLeakyLmsAvg{iLeak, iStep} = mean(cat(3, weightLeakyLms{iLeak, iStep, :}), 3);
    end
end

figure;
for iLeak = 1: nLeaks
    legendStr = cell(2 * orderAr, 1);
    for iStep = 1: nSteps
        subplot(nLeaks, nSteps, (iLeak - 1) * nSteps + iStep);
        for iOrder = 1: orderAr
            plot(weightLeakyLmsAvg{iLeak, iStep}(iOrder, :), 'LineWidth', 2);
            hold on;
            legendStr{2 * iOrder - 1} = sprintf('Est. a_%d', iOrder);
            plot([0 nSamples], [coefAr(iOrder) coefAr(iOrder)], '--', 'LineWidth', 2);
            hold on;
            legendStr{2 * iOrder} = sprintf('a_%d', iOrder);
        end
        hold off;
        grid on; grid minor;
        legend(legendStr, 'location', 'bestoutside');
        title(sprintf('Steady state values of coefficients for \\mu = %.2f and \\gamma = %.1f', step(iStep), leak(iLeak)));
        xlabel('Number of iterations (sample)');
        ylabel('Average weights');
        ylim([0 1]);
    end
end