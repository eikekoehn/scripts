%% BEC: Si/C ratio of diatom Si uptake
%
% Si/C ratio depends on:
%   surrounding Fe levels
%   surrounding Si levels
%   half-saturation constant of Si 
%   gQsi_min 
%   gQsi_max
%   cksi
%
% generally: 
%   the lower Fe, the higher the Si:C ratio
%   the lower Si, the lower the Si:C ratio
%   ratio varies between gQsi_min and gQsi_max
%
%   baseline ratio is gQsi0=0.137, which corresponds to 16:117, i.e. N:C in
%   Redfield, i.e. Si:N=1
%
%   gQsi_max = 0.685 corresponds to Si:N=5 (SO diatoms often have this or
%   even higher!)
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!

close all; clear all; clc;

% choose parameters (each combination of half-sat constants will be plotted)
kSi      = [0.8,1];  % mmol m-3
kFe      = [0.12*10^-3,0.12*10^-3]; % mmol m-3

gQsi0     = 0.137;
gQsi_max  = 0.685;
gQsi_min  = 0.0457;

% choose Si & fe levels
n = 100;
si = linspace(0,15,n); % mmol m-3
fe = linspace(0,3,n); fe = fe.*10^-3;  % mmol m-3

%%%%%
%% nothing to be changed below this line
%%%%%

cksi     = 5;

SiCratio = NaN(numel(fe),numel(si),numel(kSi));
for pp=1:numel(kSi)
    
    SiCratio(:,:,pp) = gQsi0*ones(numel(fe),numel(si));
    
    for i=1:numel(fe)
        for j=1:numel(si)
            if fe(i)<cksi*kFe(pp) && fe(i)>0 && si(j)>cksi*kSi(pp)
                
                dd = SiCratio(i,j,pp)*cksi*kFe(pp)/fe(i);
                
                if dd<gQsi_max
                    SiCratio(i,j,pp)=dd;
                else
                    SiCratio(i,j,pp)=gQsi_max;
                end
                clear dd
            end
            if fe(i)==0
                SiCratio(i,j,pp) = gQsi_max;
            end
            if si(j)<cksi*kSi(pp)
                dd = SiCratio(i,j,pp)*si(j)/(cksi*kSi(pp));
                if dd<gQsi_min
                    SiCratio(i,j,pp)=gQsi_min;
                else
                    SiCratio(i,j,pp)=dd;
                end
                clear dd
            end
        end
    end
    clear i j
    
end
clear pp



%%%%
%% plot
%%%%

width  = 22; 
height = 17;
vis    = 'on';
fs     = 14;
lw     = 3;

for  i=1:numel(kFe)
    if ~exist('siC_string','var')
        siC_string = {strcat(['kFe=',num2str(kFe(i)),...
            ', kSi=',num2str(kSi(i))])};
    else
        siC_string = [siC_string,strcat(['kFe=',num2str(kFe(i)),...
            ', kSi=',num2str(kSi(i))])];
    end
end
clear i


for ss = 1:numel(kFe)
    
    fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
    pcolor(squeeze(SiCratio(:,:,ss)))
    shading flat
    hh=colorbar;
    pos=get(hh,'Position');
    caxis([gQsi_min gQsi_max])
    text(numel(fe)+7,1,'gQsi min','FontSize',fs-2)
    text(numel(fe)+7,numel(si),'gQsi max','FontSize',fs-2)
    title(siC_string{ss},'FontSize',fs)
    ylabel('Fe [\mumol m^{-3}]','FontSize',fs)
    xlabel('Si [mmol m^{-3}]','FontSize',fs)
    set(gca,'YTick',1:10:numel(fe),'YTickLabel',round(10.*1000.*fe(1:10:end))./10)
    set(gca,'XTick',1:10:numel(si),'XTickLabel',round(10.*si(1:10:end))./10)
    set(gca,'FontSize',fs)
    
end
clear ss


