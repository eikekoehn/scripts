%% BEC: denitrification vs aerobis respiration (water column)
%
% partitioning depends on:
%   surrounding O2 levels
%   parm_o2_min
%   parm_o2_min_delta
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!
%

%ccc

% choose constants
parm_o2_min = 4;
parm_o2_min_delta = 2;

% choose O2 levels
O2 = 0:0.1:10;

%%%%%
%% nothing to be changed below this line
%%%%%

% DENITRIF
WORK1a = ((parm_o2_min + parm_o2_min_delta) - O2)./parm_o2_min_delta;
WORK1a = min(max(WORK1a,0),1);

% O2_CONS
WORK1b = (O2 - parm_o2_min) ./ parm_o2_min_delta;
WORK1b = min(max(WORK1b,0),1);


%%%%
%% plot
%%%%

colors = {'k','b','r','g','m','c','y'};
width  = 22; 
height = 17;
vis    = 'on';
fs     = 14;
lw     = 3;

fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
plot(O2,WORK1a,'b','LineWidth',lw)
hold on; box on;
plot(O2,WORK1b,'r','LineWidth',lw)
legend('Denitrif','Aerobic respiration')
xlabel('O2 [mmol m^{-3}]','FontSize',fs)
ylabel('Partitioning of total respiration','FontSize',fs)
set(gca,'FontSize',fs)




