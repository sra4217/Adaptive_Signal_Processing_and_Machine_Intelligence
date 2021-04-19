for rank=1:length(svdX)
    [U,S,V] = svd(Xnoise);
    
    U = U(:,1:rank);
    S = S(1:rank,1:rank);
    V = V(:,1:rank);

    Xrecon = U*S*V';

    e1(rank) = sum(sum((X-Xrecon).^2));
    e2(rank) = sum(sum((Xrecon-Xnoise).^2));
end

figure
plot(e1, 'linewidth', 2)
hold on
plot(e2, 'linewidth', 2)
xlabel('Rank')
ylabel('Error')
legend('X - Xrecon','Xrecon - Xnoise')
title('Squared Error vs reconstruction rank')
