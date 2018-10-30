%% BEC: Si/C ratio of diatom Si uptake
%
% Dale 2015
%
% flux of iron from sediments depends on:
%   bottom water O2 levels
%   POC flux to the sediment
%
% generally: 
%   the higher O2, the lower the Fe flux
%   the higher the POC flux, the higher the Fe flux
% Note that what I call "POC flux" here is "carbon oxidation rate" in the
% original paper. Here, the "POC flux" is corrected for what is lost to the
% sediment in BEC, and I only consider the part of the "POC flux to the
% sediment" that is actually respired. 
%
% Note that the data to fit the equation used are not really global...
% For example, for the SO, the benthic fluxes are quite high...
% Double-check iron inventory in your domain (especially scavenging rates as a sink)!!
%
% 
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!

close all; clear all; clc;

% choose maximum flux of iron (in Dale et al., gamma is equal to 170 mumol
% m-2 d-1)
gamma = 170;  % umol m-2 d-1

% choose levels of O2 and POC flux to look at
n = 100;
c_ox  = linspace(0,15,n);   % mmol m-2 d-1
o2_bw = linspace(0,200,n);  % uM


%%%%%
%% nothing to be changed below this line
%%%%%

dfe = NaN(n,n);
for i=1:n
    for j=1:n
        dfe(i,j) = gamma .* tanh(c_ox(i)./o2_bw(j));
    end
end
clear i j

%%%%
%% plot
%%%%

width  = 22; 
height = 17;
vis    = 'on';
fs     = 14;
lw     = 3;

fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
pcolor(dfe)
shading flat
hh=colorbar;
set(get(hh,'ylabel'),'string','\mumol m^{-2} d^{-1}','FontSize', fs);
caxis([0 gamma])
title(strcat(['Benthic Iron Flux, gamma=',num2str(gamma),' \mumol m^{-2} d^{-1}']),'FontSize',15,'FontWeight','bold')
xlabel('bottom water oxygen [\muM]','FontSize',14)
ylabel('carbon oxidation rate [mmol m^{-2} day^{-1}]','FontSize',14)
set(gca,'YTick',1:10:numel(c_ox),'YTickLabel',round(10.*c_ox(1:10:end))./10)
set(gca,'XTick',1:10:numel(o2_bw),'XTickLabel',round(10.*o2_bw(1:10:end))./10)
set(gca,'FontSize',fs)


