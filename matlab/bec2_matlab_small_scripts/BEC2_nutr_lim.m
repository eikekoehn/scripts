%% BEC: nutrient limitation of phytoplankton growth
%
% Michaelis Menten
%
% nutrient limitation (n.d.) depends on:
%   nutrient concentrations
%   nutrient half-saturation constants 
%
% NOTE: this script just visualizes the function used to calculate the
%   limitation. In BEC, this is done for every nutrient, and only the most-limiting nutrient
%   is used for nutrient uptake.
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!

close all; clear all; clc

% choose half-saturation constant (up to seven currently possible)
k  = [0.1,0.5,0.8];

% choose nutrient range (note that Fe is typically in the order of 10^-3 mmol
% m-3, while all others are in mmol m-3)
nutrient = 0:0.01:10;

%%%%%
%% nothing to be changed below this line
%%%%%

nutr_lim = NaN(length(nutrient),length(k));
for i=1:length(k)
    nutr_lim(:,i) = nutrient./(k(i)+nutrient);
end
clear i


%%%%
%% plot
%%%%

colors = {'k','b','r','g','m','c','y'};
width  = 22;
height = 17;
vis    = 'on';
fs     = 14;
lw     = 3;

for  i=1:numel(k)
    if ~exist('k_string','var')
        k_string = {strcat(['k=',num2str(k(i)),' mmol m^{-3}'])};
    else
        k_string = [k_string,strcat(['k=',num2str(k(i)),' mmol m^{-3}'])];
    end
end
clear i

fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
for i=1:length(k)
    plot(nutrient,nutr_lim(:,i),colors{i},'LineWidth',lw)
    hold on; box on;
end
legend(k_string,'Location','SouthEast')
xlabel('Nutrient Concentration [mmol m^{-3}]','FontSize',fs)
ylabel('Nutrient limitation [n.d.]','FontSize',fs)
set(gca,'FontSize',fs)




