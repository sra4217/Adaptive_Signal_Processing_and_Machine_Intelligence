load('./Data/EEG_Data_Assignment1.mat')

[PSDPOz,f4] = periodogram(POz,hamming(length(POz)),5*fs,fs);
figure
plot(f4,10*log10(PSDPOz), 'linewidth', 2)
xlabel('Frequency')
ylabel('Magnitude')
xlim([0 70])
title('EEG data periodogram')

POz1 = reshape(POz-mean(POz), fs, []);
POz5 = reshape(POz-mean(POz), fs*5, []);
POz10 = reshape(POz-mean(POz), fs*10, []);
PSDPOz1 = mean(periodogram(POz1,hamming(length(POz1)),5*fs,fs),2);
PSDPOz5 = mean(periodogram(POz5,hamming(length(POz5)),5*fs,fs),2);
PSDPOz10 = mean(periodogram(POz10,hamming(length(POz10)),5*fs,fs),2);

figure
plot(f4,10*log10(PSDPOz1), 'linewidth', 2);
hold on
plot(f4,10*log10(PSDPOz5), 'linewidth', 2);
plot(f4,10*log10(PSDPOz10), 'linewidth', 2);
xlim([0 70]);
xlabel('Frequency')
ylabel('Magnitude')
legend('1s window','5s window','10s window','Location','northwest')
title('Periodogram for windowed EEG data')