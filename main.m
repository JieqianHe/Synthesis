N = 128;
X = zeros(N, 1);
X(43) = 1;
%S = [1:N];
S = [1/2:1/2:log2(N)];
S = floor(2.^S);
xi = 41/N * pi;
k = 14;
nrow = 3;
[X_bar, SX_bar, Sx] = gaborwave_synth_greedy_1d(X, S, xi);

%visualize X and X_bar
figure();
subplot(nrow,(k+1)/nrow,1);
plot(X,'LineWidth',2)
%leg(1) = 'orig'
title('original');
set(gca,'FontSize',18);

for ind = 1:k
subplot(nrow,(k+1)/nrow,ind+1);
plot(X_bar(:,ind),'LineWidth',2);
title(['X-bar-' num2str(ind) ]);
set(gca,'FontSize',18);
end
%{
saveas(gca,['dirac1_512_test4.jpg']);
save('SX_bar_dirac1_512_test3.mat','SX_bar');
save('X_bar_dirac1_512_test3.mat','X_bar');
%}
figure();
subplot(1,2,1);
%subplot(3,4,1);
plot(X,'LineWidth',2)
%leg(1) = 'orig'
title('original');
set(gca,'FontSize',18);
subplot(1,2,2);
%subplot(3,4,ind+1);
plot(X_bar(:,end),'LineWidth',2);
title(['X-bar-end' ]);
set(gca,'FontSize',18);
saveas(gca,['dirac3_128_test15.png']);