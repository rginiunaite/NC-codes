
sim1 = 'aqp_M1_prop_box_plot.csv';
M1 = csvread(sim1);
M1 = M1(:,1);

sim2 = 'aqp_M2_prop_box_plot.csv';
M2 = csvread(sim2);
M2 = M2(:,1);

sim3 = 'aqp_M3_prop_box_plot.csv';
M3 = csvread(sim3);
M3 = M3(:,1);


figure
h = boxplot([M1,M2,M3],'Labels',{'Model 1','Model 2','Model 3'},'Whisker',1);
set(h,{'linew'},{3})

set(gca,'linew',3)

ylim([0 1])
yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
yticklabels({'0.0', '0.2', '0.4', '0.6', '0.8', '1.0'})
ylabel('Fraction of follower cells not in chains')

set(gca,'FontSize',36)
ax = gca;



