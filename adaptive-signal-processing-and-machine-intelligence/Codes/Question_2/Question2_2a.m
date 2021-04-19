N = 1000; 
realisations = 100; 

b = [1 0.9];
a = 1;
order = length(b);

rho = 0.005;
mu = [0.01, 0.1];

A = zeros(5, order, N, realisations);
error = zeros(5, N, realisations);
for i=1:realisations
    eta = normrnd(0, sqrt(0.5), 1, N);
    x = filter(b,a,eta);
    
    [xhat, error(1, :, i), A(1, :, :, i)] = lmsG(x, eta, mu(1), order);
    [xhat, error(2, :, i), A(2, :, :, i)] = lmsG(x, eta, mu(2), order);
    [xhat, error(3, :, i), A(3, :, :, i)] = gassLMS(x, eta, 0, order, rho, 'ben');
    [xhat, error(4, :, i), A(4, :, :, i)] = gassLMS(x, eta, 0, order, rho, 'ang');
    [xhat, error(5, :, i), A(5, :, :, i)] = gassLMS(x, eta, 0, order, rho, 'mat');
end

error = mean(error,3);
A = mean(A,4);

realA = b(2)*ones(1000,1);
wE1 = realA - squeeze(A(1,2,:));
wE2 = realA - squeeze(A(2,2,:));
wE3 = realA - squeeze(A(3,2,:));
wE4 = realA - squeeze(A(4,2,:));
wE5 = realA - squeeze(A(5,2,:));

figure
plot(wE1,'linewidth',2)
hold on
plot(wE2,'linewidth',2)
plot(wE3,'linewidth',2)
plot(wE4,'linewidth',2)
plot(wE5,'linewidth',2)
legend('LMS 0.01','LMS 0.1','ben','ang','mat')
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Weight Error Curves', 'Fontsize', 35)

figure
plot(pow2db(wE1.^2),'linewidth',2)
hold on
plot(pow2db(wE2.^2),'linewidth',2)
plot(pow2db(wE3.^2),'linewidth',2)
plot(pow2db(wE4.^2),'linewidth',2)
plot(pow2db(wE5.^2),'linewidth',2)
legend('LMS 0.01','LMS 0.1','ben','ang','mat')
xlabel('n')
ylabel('Error')
set(gca, 'Fontsize', 22)
title('Squared Error Curves', 'Fontsize', 35)
