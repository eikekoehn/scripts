%% BEC: light limitation of phytoplankton growth
%
% Geider model
%
% light limitation (n.d.) depends on:
%   PAR
%   nutrient limitation (i.e. surrounding nutrient concentrations and
%       half-sat constants)
%   temperature limitation (i.e. surrounding temp + Q10)
%   Chl:C ratio of phytoplankton type
%   alphaPI (initial slope of Photosynthesis-Irradiance curve)
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!


close all; clear all; clc

% choose alphaPIs to visualize (up to seven currently possible)
alphaPI = [0.3,0.45,0.6,0.2];

% choose parameters of PFT:
PCref   = 4.6;   % max growth rate
Q10     = 1.55;  % sensitivity in temperature limitation
thetaC  = 0.2;   % Chl/C ratio

% set surrounding environmental conditions
nut_lim = 0.1;  % e.g. 0.1=severely nutrient limited, 1=no nutr limitation
temp    = 5;    % surorunding water temperature

% choose PAR range to look at
PAR_lay = 0:1:50;

%%%%%
%% nothing to be changed below this line
%%%%%

c1      = 1;
c10     = 10;
Tref    = 30;
epsTinv = 3.17e-08;

% temp lim
Tfunc   = Q10^(((temp + 273.15)-(Tref + 273.15))/c10);

% temp & nutr limited growth
PCmax   = PCref .* nut_lim .* Tfunc;  % growth rate corrected for nutrient levels and temperature

light_lim = NaN(length(PAR_lay),length(alphaPI));
for i=1:length(alphaPI)
    light_lim(:,i) = c1 - exp((-c1.*alphaPI(i).*thetaC.*PAR_lay)./(PCmax+epsTinv));
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

for  i=1:numel(alphaPI)
    if ~exist('alphaPI_string','var')
        alphaPI_string = {strcat(['\alpha=',num2str(alphaPI(i))])};
    else
        alphaPI_string = [alphaPI_string,strcat(['\alpha=',num2str(alphaPI(i))])];
    end
end
clear i


fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
for i=1:length(alphaPI)
    plot(PAR_lay,light_lim(:,i),colors{i},'LineWidth',lw)
    hold on;
end
grid on;
legend(alphaPI_string,'Location','SouthEast')
xlabel('PAR [W m^{-2}]','FontSize',fs)
ylabel('light limitation [n.d.]','FontSize',fs)
title(strcat(['nutlim=',num2str(nut_lim),', temp=',num2str(temp),', PCref=',num2str(PCref),...
    ',Q10=',num2str(Q10),',thetaC=',num2str(thetaC)]),'FontSize',12)
set(gca,'FontSize',fs)


