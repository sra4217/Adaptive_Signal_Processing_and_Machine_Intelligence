fSample = 1500;
nSamples = 1500;
variance = 0.05;
nSegs = 3;
lengthSeg = 500;
freqFunc = @(n) ((1 <= n) & (n <= 500)) .* 100 + ((501 <= n) & (n <= 1000)) .* (100 + (n - 500) / 2) + ((1001 <= n) & (n <= 1500)) .* (100 + ((n - 1000) / 25) .^ 2);
freqSeq = freqFunc(1: nSamples);
phaseSeq = cumsum(freqSeq);
orderAr = 1;

coefArEstSeg = cell(nSegs, 1);
hArSeg = cell(nSegs, 1);
fArSeg = cell(nSegs, 1);
psdArSeg = cell(nSegs, 1);
fmSignal = exp(1i * 2 * pi / fSample * phaseSeq) + sqrt(variance / 2) * (randn(1, nSamples) + 1i * randn(1, nSamples));
coefArEst = aryule(fmSignal, orderAr);
[hAr, fAr] = freqz(1, coefArEst, nSamples, fSample);
psdAr = abs(hAr) .^ 2;
for iSeg = 1: nSegs
    coefArEstSeg{iSeg} = aryule(fmSignal((iSeg - 1) * lengthSeg + 1: iSeg * lengthSeg), orderAr);
    [hArSeg{iSeg}, fArSeg{iSeg}] = freqz(1, coefArEstSeg{iSeg}, nSamples / nSegs, fSample);
    psdArSeg{iSeg} = abs(hArSeg{iSeg}) .^ 2;
end

figure;
subplot(2, 1, 1);
plot(freqSeq);
grid on; grid minor;
legend('Frequency');
title('Frequency of FM signal');
xlabel('Time (sample)');
ylabel('Frequency (Hz)');
subplot(2, 1, 2);
plot(angle(exp(1i * 2 * pi / fSample * phaseSeq)));
grid on; grid minor;
legend('Phase');
title('Phase of FM signal');
xlabel('Time (sample)');
ylabel('Angle (rad)');
figure;
plot(fAr, pow2db(psdAr), 'LineWidth', 2);
grid on; grid minor;
legend(sprintf('Aryule-AR (%d)', orderAr));
title('Overall estimation of FM signal by AR model');
xlabel('Frequency (Hz)');
ylabel('PSD (dB)');
figure;
for iSeg = 1: nSegs
    subplot(nSegs, 1, iSeg);
    plot(fArSeg{iSeg}, pow2db(psdArSeg{iSeg}), 'LineWidth', 2);
    grid on; grid minor;
    legend(sprintf('Segment %d', iSeg));
    title(sprintf('Individual estimation of FM signal by Aryule-AR (%d) model', orderAr));
    xlabel('Frequency (Hz)');
    ylabel('PSD (dB)');
end