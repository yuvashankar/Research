close all;
%Execution Time

figure
hold on
title ('Speed Up Versus Threads');
xlabel('Threads');
ylabel('Speed Up');
plot (bic(2:end,1), bic(2:end,3), '-*k');

figure
hold on
title('Efficiency Versus Threads');
xlabel('Threads');
ylabel('Efficiency');
ylim([0 1]);
plot (bic(2:end,1), bic(2:end,4), '-*k');