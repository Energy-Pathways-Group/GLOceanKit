XX = 0.5:0.1:12;
load SEM1.mat
a2a=figure; plot(XX,(sum(SEM>0.5)+1)/48,'b','linewidth',2);
load SEM3.mat
hold on; plot(XX,(sum(SEM>0.5))/24,'r','linewidth',2);
load SEM4.mat
hold on; plot(XX,(sum(SEM>0.5))/24,'g','linewidth',2);
xlabel('1/\sigma (days)'); ylabel('window length (days)');
ylim([2/48 6]);
exportfig(a2a, 'figures/append2a.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);

load SEM1.mat
a2b=figure; plot(XX,(sum(SEM>0.5)+1)/48,'b','linewidth',2);
load SEM5.mat
hold on; plot(XX,(sum(SEM>0.5))/24,'r','linewidth',2);
load SEM6.mat
hold on; plot(XX,(sum(SEM>0.5))/24,'g','linewidth',2);
xlabel('1/\sigma (days)'); ylabel('window length (days)');
ylim([2/48 6]);
exportfig(a2b, 'figures/append2b.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);