clear all
%%%% average simulations

simulations1 = 'aqp_M1.csv';
M1 = csvread(simulations1); % 
M1 = M1(:,1)/10;

simulations2 = 'aqp_M2.csv';
M2 = csvread(simulations2);
M2 = M2(:,1)/10;

simulations3 = 'aqp_M3.csv';
M3 = csvread(simulations3); 
M3 = M3(:,1)/10;
 
simulations4 = 'aqp_M4.csv';
M4 = csvread(simulations4); 
M4 = M4(:,1)/10;

simulations5 = 'aqp_M5.csv';
M5 = csvread(simulations5); 
M5 = M5(:,1)/10;

simulations6 = 'aqp_M6.csv';
M6 = csvread(simulations6); 
M6 = M6(:,1)/10;

simulations7 = 'aqp_M7.csv';
M7 = csvread(simulations7); 
M7 = M7(:,1)/10;

simulations8 = 'aqp_M8.csv';
M8 = csvread(simulations8); 
M8 = M8(:,1)/10;

simulations9 = 'aqp_M9.csv';
M9 = csvread(simulations9); 
M9 = M9(:,1)/10;

x = [0:50:1000];


figure 


options = fitoptions('Method','Smooth','SmoothingParam',0.0000011);

%%% Plot data
[f1,gof,out] = fit(x',M1,'smooth',options);
[f2,gof,out] = fit(x',M2,'smooth',options);
[f3,gof,out] = fit(x',M3,'smooth',options);
[f4,gof,out] = fit(x',M4,'smooth',options);
[f5,gof,out] = fit(x',M5,'smooth',options);
[f6,gof,out] = fit(x',M6,'smooth',options);
[f7,gof,out] = fit(x',M7,'smooth',options);
[f8,gof,out] = fit(x',M8,'smooth',options);
[f9,gof,out] = fit(x',M9,'smooth',options);
%plot(f1, x',M1) % with data points

 h1 = plot(f1,'blue')
hold on
h2 = plot(f2,'red')
h3 = plot(f3,'green')

%h4 = plot(f4,'black')
%h5 = plot(f5,'green')
% h6 = plot(f6,'magenta')
% h7 = plot(f7,'cyan')
% h8 = plot(f8,'black')
% h9 = plot(f9,'magenta')


legend ('Model 1','Model 2','Model 3','M4','M5','M6','M7')
legend ('Model 1','Model 2','Model 3')


% hold off
% figure
% h8 = plot(f8, 'red')
% hold on
% h9 = plot(f9,'blue')

ylim([0 20])
h1.LineWidth =3;
h2.LineWidth =3;
h3.LineWidth =3;
 h4.LineWidth =3;
 h5.LineWidth =3;
 h6.LineWidth =3;
 h7.LineWidth =3;
 h8.LineWidth =3;
 h9.LineWidth = 3;



xlabel('Distance from the neural tube, \mu m','FontSize',14)
set(gca,'linewidth',2)
ylabel('Number of cells','FontSize',14)
set(gca,'FontSize',36)
 ax = gca;




% 
% [f2,gof,out] = fit(x',M2,'smooth',options);
% 
% 
% hold on
% plot(f1, x',M_AQup,'blue')

 