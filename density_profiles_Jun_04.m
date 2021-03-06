clear all
%%%% average simulations

simulations1 = 'AQP-1-GOOD/aqp_M1_Jul5.csv';

M1 = csvread(simulations1); % 
M1 = M1(:,1)/10;

simulations2 = 'AQP-1-GOOD/aqp_M2_Jul5.csv';
M2 = csvread(simulations2);
M2 = M2(:,1)/10;
% 
 simulations3 = 'AQP-1-GOOD/aqp_M3_Jul5.csv';
M3 = csvread(simulations3); 
M3 = M3(:,1)/10;
 
simulations4 = 'AQP-1-GOOD/aqp_M4_Jul5.csv';
M4 = csvread(simulations4); 
M4 = M4(:,1)/10;

simulations5 = 'AQP-1-GOOD/aqp_M5_Jul5.csv';
M5 = csvread(simulations5); 
M5 = M5(:,1)/10;

simulations6 = 'AQP-1-GOOD/aqp_M6_Jul5.csv';
M6 = csvread(simulations6); 
M6 = M6(:,1)/10;
% 
simulations7 ='AQP-1-GOOD/aqp_M7_Jul5.csv';
M7 = csvread(simulations7); 
M7 = M7(:,1)/10;
% 
simulations8 = 'AQP-1-GOOD/aqp_M8_Jul5.csv';
M8 = csvread(simulations8); 
M8 = M8(:,1)/10;

simulations9 = 'AQP-1-GOOD/aqp_M9_Jul5.csv';
M9 = csvread(simulations9); 
M9 = M9(:,1)/10;

 simulations10 = 'AQP-1-GOOD/aqp_M10_Jul5.csv';

 M10 = csvread(simulations10);
 M10 = M10(:,1)/10;

x = [0:55:1100];


 simulations11 = 'AQP-1-GOOD/aqp_M11_Jul5.csv';

 M11= csvread(simulations11);
 M11 = M11(:,1)/10;

x = [0:55:1100];

figure 


options = fitoptions('Method','Smooth','SmoothingParam',0.0000011);

%%% Plot data
%[f1,gof,out] = fit(x',M1,'smooth',options);
%[f2,gof,out] = fit(x',M2,'smooth',options);
%[f3,gof,out] = fit(x',M3,'smooth',options);
% [f4,gof,out] = fit(x',M4,'smooth',options);
% [f5,gof,out] = fit(x',M5,'smooth',options);
% [f6,gof,out] = fit(x',M6,'smooth',options);
% [f7,gof,out] = fit(x',M7,'smooth',options);
% [f8,gof,out] = fit(x',M8,'smooth',options);
% [f9,gof,out] = fit(x',M9,'smooth',options);
%plot(f1, x',M1) % with data points

%h1 = plot(f1,'blue')
hold on
%h2 = plot(f2,'red')
%h3 = plot(f3,'green')

%h4 = plot(f4,'black')
%h5 = plot(f5,'green')
% h6 = plot(f6,'magenta')
% h7 = plot(f7,'cyan')
% h8 = plot(f8,'black')
% h9 = plot(f9,'magenta')
 
 
  h1 = plot(x,M1,'--k')
  h2 = plot(x,M2,':','Color',[0.6350 0.0780 0.1840])
%  h3 = plot(x,M3,'Color',[0.9290 0.6940 0.1250])
% h4 = plot(x,M4,'Color',[0.4660 0.6740 0.1880])
 %h5 = plot(x,M5,':','Color',[0.3010 0.7450 0.9330])
 %h6 = plot(x,M6,'-.','Color', '[0.6350 0.0780 0.1840]')
% h7 = plot(x,M7,'--','Color',[0 0.4470 0.7410])
%  h8 = plot(x,M8,'-.','Color',[0.9290 0.6940 0.1250])
%  h9 = plot(x,M9,'Color',[0.4940 0.1840 0.5560])
  h10 = plot(x,M10,'-.','Color',[0.3010 0.7450 0.9330])
  h11 = plot(x,M11,'Color',[0.4660 0.6740 0.1880])

 
 %h1 = stairs(x,M1)
 %h2 = stairs(x,M2)
 %h3 = stairs(x,M3)
 %h4 = stairs(x,M4)
 %h5 = stairs(x,M5)
 %h6 = stairs(x,M6)
 %h7 = stairs(x,M7)
 %h8 = stairs(x,M8)
 %h9 = stairs(x,M9)


%legend ('Model 1','Model 2','Model 3')%,'M4','M5','M6','M7')
%legend('Model 4','Model 5','Model 6','Model 7')
%legend('Model 2', 'Model 7','Model 8','Model 9')
legend ('Model 1','Model 2','Model 10','Model 11')


% hold off
% figure
% h8 = plot(f8, 'red')
% hold on
% h9 = plot(f9,'blue')

%ylim([0 20])
h1.LineWidth =3;
h2.LineWidth =3;
h3.LineWidth =3;

h4.LineWidth =3;
 h5.LineWidth =3;
 h6.LineWidth =3;
 h7.LineWidth =3;
 h8.LineWidth =3;
 h9.LineWidth = 3;
 h10.LineWidth = 3;
 h11.LineWidth = 3;

xlim([0 1100])


xlabel(['Distance from the neural tube, ',char(181),'m'],'FontSize',14)
set(gca,'linewidth',2)
ylabel(['Number of cells (per 55 ',char(181),'m)'],'FontSize',14)
set(gca,'FontSize',36)
ax = gca;
box on
grid on


% 
% [f2,gof,out] = fit(x',M2,'smooth',options);
% 
% 
% hold on
% plot(f1, x',M_AQup,'blue')

 