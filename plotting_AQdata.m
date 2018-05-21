clear all
%%%% average simulations

simulations = 'aqp_down.csv';
M_AQdown = csvread(simulations);
M_AQdown = M_AQdown(:,1)/10;

simulations2 = 'control.csv';
M_control = csvread(simulations2); % 
M_control = M_control(:,1)/10;

simulations1 = 'aqp_up.csv';
M_AQup = csvread(simulations1); 
M_AQup = M_AQup(:,1)/10;
 
 simulations3 = 'aqp_up_slower.csv';
 M_AQup_dir = csvread(simulations3); % 
 M_AQup_dir = M_AQup_dir(:,1)/10;
 
Matrix_comb(:,1) = M_control;
Matrix_comb(:,2) = M_control;

Matrix_comb2(:,1) = M_AQdown;
Matrix_comb2(:,2) = M_AQdown;




Matrix_comb3(:,1) = M_AQup;
Matrix_comb3(:,2) = M_AQup;

 Matrix_comb4(:,1) = M_AQup_dir;
 Matrix_comb4(:,2) = M_AQup_dir;
 
X = [1:2];

Y = [0:20];

figure
% plot underexpression
s = surf(X,Y,Matrix_comb)
colorbar
shading interp

hold on 
% plot control
X2 = [2,3];
s = surf (X2,Y, Matrix_comb2)
shading interp

% plot overexpression
X3 = [3,4];
s = surf (X3,Y, Matrix_comb3)
shading interp

% directional overexpression

 X4 = [4,5];
 s= surf(X4, Y, Matrix_comb4)
 shading interp


yticks([0 3 6 9 12 15 18 21])
yticklabels({'0','175','350','525','700', '875','950','1175'})

ylabel('\mum', 'FontSize',12);

set(gca,'xtick',[])
set(gca,'xticklabel',[])

set(gca,'FontSize',24)

t = annotation('textbox');
sz = t.FontSize;
t.FontSize = 24;
