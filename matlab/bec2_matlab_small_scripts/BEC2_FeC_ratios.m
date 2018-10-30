%% BEC: Fe/C ratio of phytoplankton fe uptake
%
% Fe/C ratio depends on:
%   surrounding Fe levels
%   half-saturation constant of iron of respective phytoplankton
%   gQfe0 (it's called "initital ratio" in code, I think it's more a
%   maximum Fe/C ratio at high dFe)
%   gQfe_min (minimum Fe/C ratio when dFe is low)
%   cks (used for onset of transition from gQfe0 to gQfe_min)
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!

%ccc

% choose parameters (up to seven combinations currently possible)
kFe      = [0.3*10^-3,0.18*10^-3,0.12*10^-3];  % mmol m-3
gQfe0    = [20*10^-6,20*10^-6,20*10^-6];  
gQfe_min = [3*10^-6,3*10^-6,3*10^-6];
cks      = [9,9,9];

% choose iron levels
fe = 0:0.005:3; fe = fe.*10^-3;  % mmol m-3


%%%%%
%% nothing to be changed below this line
%%%%%

FeCratio = NaN(numel(fe),numel(kFe));
for pp=1:numel(kFe)
    
    FeCratio(:,pp) = gQfe0(pp)*ones(numel(fe),1);
    
    for i=1:numel(fe)
        if fe(i)<cks(pp)*kFe(pp)
            dd = FeCratio(i,pp)*fe(i)/(cks(pp)*kFe(pp));
            if dd>gQfe_min(pp)
                FeCratio(i,pp)=dd;
            else
                FeCratio(i,pp)=gQfe_min(pp);
            end
            clear dd
        end
    end
    clear i
    
end
clear pp



%%%%
%% plot
%%%%

colors = {'k','b','r','g','m','c','y'};
width  = 22; 
height = 17;
vis    = 'on';
fs     = 14;
lw     = 3;

for  i=1:numel(kFe)
    if ~exist('feC_string','var')
        feC_string = {strcat(['kFe=',num2str(kFe(i)),...
            ', gQfe0=',num2str(gQfe0(i)),', gQfe min=',num2str(gQfe_min(i)),', cks=',num2str(cks(i))])};
    else
        feC_string = [feC_string,strcat(['kFe=',num2str(kFe(i)),...
            ', gQfe0=',num2str(gQfe0(i)),', gQfe min=',num2str(gQfe_min(i)),', cks=',num2str(cks(i))])];
    end
end
clear i


fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
for i=1:length(kFe)
    plot(fe,FeCratio(:,i),colors{i},'LineWidth',lw)
    hold on;
end
grid on;
legend(feC_string,'Location','SouthEast')
xlabel('Fe [mmol m^{-3}]','FontSize',fs)
ylabel('Fe/C ratio [mmol Fe (mmol C)^{-1}]','FontSize',fs)
set(gca,'FontSize',fs)







