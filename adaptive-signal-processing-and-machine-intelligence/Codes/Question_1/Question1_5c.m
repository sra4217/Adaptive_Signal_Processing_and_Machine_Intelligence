ecg = load('./data/ECG_Data_Assignment1.mat');
fSample = ecg.fsRRI;
rri = {detrend(ecg.xRRI1 - mean(ecg.xRRI1)) detrend(ecg.xRRI2 - mean(ecg.xRRI2)) detrend(ecg.xRRI3 - mean(ecg.xRRI3))};
nRris = length(rri);
label = ["normal", "fast", "slow"];
nFft = 2048;
orderAr = 2: 4: 10;
nOrders = length(orderAr);

psdStd = cell(nRris, 1);
for iRri = 1: nRris
    nSamples = length(rri{iRri});
    [psdStd{iRri}, f] = periodogram(rri{iRri}, rectwin(nSamples), nFft, fSample);
end

psdAr = cell(nRris, nOrders);
varEst = zeros(nRris, nOrders);
fAr = cell(nRris, 1);
for iRri = 1: nRris
    nSamples = length(rri{iRri});
    for iOrder = 1: nOrders
        [coefArEst, varEst(iRri, iOrder)] = aryule(rri{iRri}, orderAr(iOrder));
        [hAr, fAr{iRri}] = freqz(sqrt(varEst(iRri, iOrder)), coefArEst, nSamples, fSample);
        psdAr{iRri, iOrder} = abs(hAr) .^ 2;
    end
end

legendStr = cell(1, nOrders + 1);
figure;
for iRri = 1: nRris
    subplot(nRris, 1, iRri);
    plot(f, pow2db(psdStd{iRri}), 'g', 'LineWidth', 1);
    hold on;
    legendStr{1} = 'Standard';
    for iOrder = 1: nOrders
        plot(fAr{iRri}, pow2db(psdAr{iRri, iOrder}), 'k', 'LineWidth', 1);
        hold on;
        legendStr{iOrder + 1} = sprintf('AR (%d)', orderAr(iOrder));
    end
    grid on; grid minor;
    legend(legendStr);
    title(sprintf('PSD estimate by normal periodogram and AR model for %s RRI', label(iRri)));
    xlabel('Frequency');
    ylabel('PSD');
    ylim([-80 0]);
end
