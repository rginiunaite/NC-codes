%%%% average simulations

simulations = 'simulations_domain_partition_simple_downregulation.csv';
M_AQdown = csvread(simulations);
M_AQdown = M_AQdown(:,1);

simulations1 = 'simulations_domain_partition_simple_persistent.csv';
M_AQup = csvread(simulations1); 
M_AQup = M_AQup(:,1);
 
simulations2 = 'simulations_domain_partition_simple_control.csv';
M_control = csvread(simulations2); % 
M_control = M_control(:,1);
 
Matrix_comb(:,1) = M_AQdown/20;
%Matrix_comb(:,2) = M_AQdown/20;


Matrix_comb(:,2) = M_AQup/20;
%Matrix_comb(:,4) = M_AQup/20;

Matrix_comb(:,3) = M_control/20;
%Matrix_comb(:,4) = M_control/20;
 
X = [1:3];

Y = [1:13];

s = surf(X,Y,Matrix_comb)
colorbar

%shading interp


% 
% step = (100*0.005-0.005)/length(M_lin(:,1));
% 
% threshold = 0.005:step:100*0.005-0.001;
% 
% figure
% surf(M_lin,'LineWidth',2);
% shading interp
% 
% hold on 
% 
% % plot(threshold,V_quad,'LineWidth',2);
% % 
% % plot(threshold,V_log,'LineWidth',2);
% % 
% % plot(threshold, V_lam100, 'LineWidth',2);
% % 
% % plot(threshold, V_lam500, 'LineWidth',2);
% % 
% % ylabel('furthest distance travelled, \mum', 'FontSize',12);
% % 
% % xlabel('Sensing accuracy','FontSize',12);
% % 
% % legend('fixed, linear','fixed, quadratic','fixed, log','cell-induced, \lambda = 100','cell-induced, \lambda = 500')
% % 
% % set(gca,'FontSize',24)
% 
% ylabel('Distance from the neural tube, \mum', 'FontSize',12);
% 
% xlabel('Sensing accuracy','FontSize',12);
% 
% 
% xticks([0 20 40 60 80 100])
% xticklabels({'0.05','0.5','1','1.5','2','2.5'})
% 
% yticks([0 2 4 6 8 10 12 14 16 18])
% yticklabels({'0','10','20','30','40','50','60','70','80','90'})
% 
% set(gca,'FontSize',24)