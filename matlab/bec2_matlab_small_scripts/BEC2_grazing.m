%% BEC: grazing on phytoplankton
%
% baseline (as of Sep 2018) is a Holling Type II formulation without
% sum(phyto) in the denominator -> see Nissen et al., 2018
%
%
% grazing rate (mmol C m-3 d-1) depends on:
%   zooplankton biomass
%   max. growth rate of zoopl. grazing on respective phyto
%   temperature limitation function of zoopl.
%   biomass of respective phytoplankton
%   half-sat. constant for ingestion of respective phytoplankton
%   
%   most importantly, grazing rate depends on chosen formulation!!!!
%       Holling-Type I (Ivlev)
%       Holling-Type II (Michaelis-Menten)
%       Holling-Type III
%           Linear scaling of grazing rate with relative importance of
%               respective phytopl. in community
%           Kill-the-winner (as above, but with non-linear scaling)
%       ...
%
% Note that BEC historically did not consider the total phyto biomass
% anywhere in the grazing equation (see e.g. Moore et al., 2002, Nissen et
% al., 2018), i.e. the grazer could be "full" when feeding on phyto 1, but
% would nevertheless continue grazing on the other phyto...
%
% Note that the last two do not impact the ingestion functin and can
% therefore be combined with any of the Holling Types!!
%
% Note that this script only visualizes the difference between the three
% Holling Types (w/o sum in denominator)
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!

close all; clear all; clc;

% choose growth rates of zoopl feeding on phyto X and half-saturation constant for ingestion 
% (up to seven combinations currently possible)
z_umax = [4.0,4.4];  % growth rate
z_grz  = [1.05,1.05];  % half-saturation constant for ingestion

% choose other constants
zooC    = 5;  % zoopl biomass (mmol C m-3)
temp    = 10; % surrounding temperature (for temp. limitation of zoopl growth)
Q10_zoo = 1.5; % Q10 of zoopl (see BEC2_temp_lim.m for sensitivity of formulation to this parameter)

% choose levels of phyto biomass to look at (mmol C m-3)
phyto = 0:0.1:10;


%%%%%
%% nothing to be changed below this line
%%%%%

% temp lim of zoopl
c10     = 10;
Tref    = 30;
temp_lim = Q10_zoo.^(((temp + 273.15)-(Tref + 273.15))/c10);

% grazing rate
graze   = NaN(numel(phyto),numel(z_umax),3); % Holling I,Holling II,Holling III 
for i=1:numel(z_umax)
    graze(:,i,1) = z_umax(i) .* zooC .* temp_lim .* phyto;
    graze(:,i,2) = z_umax(i) .* zooC .* temp_lim .*(phyto./(phyto + z_grz(i)));
    graze(:,i,3) = z_umax(i) .* zooC .* temp_lim .*((phyto.^2)./((phyto.^2) + (z_grz(i).^2)));
end
clear i


%%%%
%% plot
%%%%

colors = {'k','b','r','g','m','c','y'};
width  = 42; 
height = 14;
vis    = 'on';
fs     = 14;
lw     = 3;

for  i=1:numel(z_umax)
    if ~exist('z_string','var')
        z_string = {strcat(['zumax=',num2str(z_umax(i)),', zgrz=',num2str(z_grz(i))])};
    else
        z_string = [z_string,strcat(['zumax=',num2str(z_umax(i)),', zgrz=',num2str(z_grz(i))])];
    end
end
clear i

fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);

subplot(1,3,1)
for i=1:length(z_umax)
    plot(phyto,graze(:,i,1),colors{i},'LineWidth',lw)
    hold on;
end
grid on;
legend(z_string,'Location','SouthEast')
title('Holling Type I (no max.)','FontSize',fs)
xlabel('phyto biomass [mmol C m^{-3}]','FontSize',fs)
ylabel('grazing rate [mmol C m^{-3} d^{-1}]','FontSize',fs)
set(gca,'FontSize',fs)

subplot(1,3,2)
for i=1:length(z_umax)
    plot(phyto,graze(:,i,2),colors{i},'LineWidth',lw)
    hold on;
end
grid on;
legend(z_string,'Location','SouthEast')
title('Holling Type II','FontSize',fs)
xlabel('phyto biomass [mmol C m^{-3}]','FontSize',fs)
ylabel('grazing rate [mmol C m^{-3} d^{-1}]','FontSize',fs)
set(gca,'FontSize',fs)

subplot(1,3,3)
for i=1:length(z_umax)
    plot(phyto,graze(:,i,3),colors{i},'LineWidth',lw)
    hold on;
end
grid on;
legend(z_string,'Location','SouthEast')
title('Holling Type III','FontSize',fs)
xlabel('phyto biomass [mmol C m^{-3}]','FontSize',fs)
ylabel('grazing rate [mmol C m^{-3} d^{-1}]','FontSize',fs)
set(gca,'FontSize',fs)
text(max(phyto)+0.2,4,strcat(['zooC=',num2str(zooC),' mmol m^{-3}']),'FontSize',fs-2)
text(max(phyto)+0.2,3.5,strcat(['temp=',num2str(temp),'°C']),'FontSize',fs-2)

