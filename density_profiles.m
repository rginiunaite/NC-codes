clear all
%%%% average simulations

simulations = 'control.csv';
M_control = csvread(simulations); % 
M_control = M_control(:,1)/10;

simulations2 = 'aqp_down.csv';
M_AQdown = csvread(simulations2);
M_AQdown = M_AQdown(:,1)/10;

simulations1 = 'aqp_up.csv';
M_AQup = csvread(simulations1); 
M_AQup = M_AQup(:,1)/10;
 
 simulations3 = 'aqp_up_slower.csv';
 M_AQup_dir = csvread(simulations3); % 
 M_AQup_dir = M_AQup_dir(:,1)/10;
 
 x= [0:50:1000];
 a = smooth(M_AQup)
 b = smooth(a)
 c = smooth(b)
 figure
%  scatter (x,M_control)
% % figure 
%  %plot (x, smooth(smooth(M_control)))
%  %hold on
%  %plot(x, smooth(smooth(M_AQup)))
%  
%  x_in = linspace(min(x),max(x))
% fu_interp = interp1(x,M_control,x_in,'spline' )
%  hold on
%  plot(x_in,fu_interp)
%  plot(x, smooth(M_control))
 
 
options = fitoptions('Method','Smooth','SmoothingParam',0.00001);
[f_control,gof,out] = fit(x',M_control,'smooth',options);

 
plot(f_control, x',M_control)

[f_AQ_up,gof,out] = fit(x',M_AQup,'smooth',options);

 hold on
plot(f_AQ_up, x',M_AQup)

 