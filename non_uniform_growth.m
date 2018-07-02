x1 = [0,3];
x2 = [3,6];
x3 = [6,10];


y1 = 1.5*x1+3;
y2 = 7.5;
y3 = 0;

figure
plot(x1,y1,'r','LineWidth',2)

hold on 
plot(x2, [y2,y2],'r','LineWidth',2)

plot(x3,[y3,y3],'r','LineWidth',2)
plot([6,6],[0,7.5],'--r','LineWidth',2)


   xlim([0 10])
   ylim([0 10])
   
   set(gca,'xtick',0)
set(gca,'xticklabel',0)

   set(gca,'ytick',0)
set(gca,'yticklabel',0)
set(gca, 'FontSize', 18)


   
   xlabel('Distance from the neural tube, \mu m','fontsize',18)
   ylabel('Domain growth rate','fontsize',18)