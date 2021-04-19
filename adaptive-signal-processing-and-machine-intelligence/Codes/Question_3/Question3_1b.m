load('data/WS_Data/low-wind.mat');
wind(1, :) = (v_east + 1i * v_north).';
load('data/WS_Data/medium-wind.mat');
wind(2, :) = (v_east + 1i * v_north).';
load('data/WS_Data/high-wind.mat');
wind(3, :) = (v_east + 1i * v_north).';
nWinds = 3;
orderFilter = 1: 24;
nOrders = length(orderFilter);
step = [1e-1, 1e-2, 1e-3];
leak = 0;

predictionClms = cell(nWinds, nOrders);
predictionAclms = cell(nWinds, nOrders);
errorClms = cell(nWinds, nOrders);
errorAclms = cell(nWinds, nOrders);
mspeClms = zeros(nWinds, nOrders);
mspeAclms = zeros(nWinds, nOrders);
circularityCoef = zeros(nWinds, 1);
for iWind = 1: nWinds
    [circularityCoef(iWind), ~] = circularity(wind(iWind, :));
    for iOrder = 1: nOrders
        [group] = preprocessing(wind(iWind, :), orderFilter(iOrder), 1);
        [~, predictionClms{iWind, iOrder}, errorClms{iWind, iOrder}] = clms(group, wind(iWind, :), step(iWind), leak);
        [~, ~, predictionAclms{iWind, iOrder}, errorAclms{iWind, iOrder}] = aclms(group, wind(iWind, :), step(iWind), leak);
        mspeClms(iWind, iOrder) = mean(abs(errorClms{iWind, iOrder}) .^ 2);
        mspeAclms(iWind, iOrder) = mean(abs(errorAclms{iWind, iOrder}) .^ 2);
    end
end

plotStep = 5;
for iWind = 1: nWinds
    fig = figure('name', sprintf('Wind level %d', iWind));
    counter = 0;
    for iOrder = 1: plotStep: nOrders
        counter = counter + 1;
        subplot(ceil(nOrders / plotStep), 2, counter);
        scatter(real(wind(iWind, :)), imag(wind(iWind, :)), 3, 'filled', 'k');
        hold on;
        scatter(real(predictionClms{iWind, iOrder}), imag(predictionClms{iWind, iOrder}), 3, 'filled', 'b');
        legend('Measured', 'CLMS', 'location', 'bestoutside');
        title(sprintf('\\rho = %.2f and M = %d', circularityCoef(iWind), orderFilter(iOrder)));
        xlabel('Real part');
        ylabel('Imaginary part');
        counter = counter + 1;
        subplot(ceil(nOrders / plotStep), 2, counter);
        scatter(real(wind(iWind, :)), imag(wind(iWind, :)), 3, 'filled', 'k');
        hold on;
        scatter(real(predictionAclms{iWind, iOrder}), imag(predictionAclms{iWind, iOrder}), 3, 'filled', 'r');
        legend('Measured', 'ACLMS', 'location', 'bestoutside');
        title(sprintf('\\rho = %.2f and M = %d', circularityCoef(iWind), orderFilter(iOrder)));
        xlabel('Real part');
        ylabel('Imaginary part');
    end
end
figure;
for iWind = 1: nWinds
    subplot(nWinds, 1, iWind);
    plot(pow2db(mspeClms(iWind, :)), 'LineWidth', 2);
    hold on;
    plot(pow2db(mspeAclms(iWind, :)), 'LineWidth', 2);
    hold off;
    grid on; grid minor;
    legend('CLMS', 'ACLMS');
    title(sprintf('Learning curves of wind level %d', iWind));
    xlabel('Order');
    ylabel('MSPE (dB)');
end
