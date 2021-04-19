
N = 10000; 
n = 1:N;
sigma = 1;
x = randn(size(n));

x(1) = x(1);
x(2) = 2.76*x(1) + x(2);
x(3) = 2.76*x(2) - 3.81*x(1) + x(3);
x(4) = 2.76*x(3) - 3.81*x(2) + 2.65*x(1) + x(4);

for i = 5:N
 x(i) = 2.76*x(i-1) - 3.81*x(i-2) + 2.65*x(i-3) - 0.92*x(i-4) + x(i);
end

x = x(5001:end);

p=2:14;
h = zeros(length(p), N/2);
aOrig = [2.76 -3.81 2.65 -0.92];
hOrig = freqz(1 , [1 -aOrig], N/2);
for i=1:length(p)
    [a,e] = aryule(x,p(i));
    [h(i,:),~] = freqz(e^(1/2),a,N/2);
end

normf = linspace(0,1,length(hOrig));
figure
plot(normf,pow2db(abs(h(2,:)').^2),'c', 'linewidth',2)
hold on
plot(normf,pow2db(abs(h(4,:)').^2),'b', 'linewidth',2)
hold on
plot(normf,pow2db(abs(h(8,:)').^2),'g', 'linewidth',2)
hold on
plot(normf,pow2db(abs(h(11,:)').^2),'y', 'linewidth',2)
hold on
plot(normf,pow2db(abs(hOrig).^2), 'k', 'linewidth',2) 
xlim([0 0.6])
xlabel('Normalised Frequency')
ylabel('Magnitude')

legend('AR(2)','AR(4)','AR(8)','AR(11)','Original')
title('Spectrum estimation for varying orders')
