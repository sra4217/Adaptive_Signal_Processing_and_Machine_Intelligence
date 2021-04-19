N = 1000; 
realisations = 100; 

b = [1 0.9];
a = 1;
order = length(b);

rho = 0.001;
mu = 0.01;

A = zeros(2, order, N, realisations);
error = zeros(2, N, realisations);
eps = ones(1,N); 
for i=1:realisations
    eta = normrnd(0, sqrt(0.5), 1, N);
    x = filter(b,a,eta);
    
    [xhat, error(1, :, i), A(1, :, :, i)] = gassLMS(x, eta, mu, order, rho, 'ben');
    [xhat, error(2, :, i), A(2, :, :, i)] = gngd(x, eta, mu, order, eps, rho);
end

error = mean(error,3);
A = mean(A,4);

realA = b(2)*ones(1,1000);
wBenError = realA - squeeze(A(1,2,:));
wGngdError = realA - squeeze(A(2,2,:));

figure
plot(wBenError(:,1), 'b','linewidth',2)
hold on
plot(wGngdError(:,1), 'k','linewidth',2)
xlim([0 200])
legend('Ben','GNGD')
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Squared Error Curves rho=0.001 mu=0.01', 'Fontsize', 35)
grid on
grid minor