%% BEC: temperature limitation of phytoplankton growth
%
% Q10 formulation
%
% temperature limitation (n.d.) depends on:
%   temperature
%   Q10 (temperature sensitivity)
%
% CN (Sep 2018): please correct any bug you find and update the files on
% the Wiki as parametrizations change!

close all; clear all; clc

% choose Q10s to visualize (up to seven currently possible)
Q10 = [1.2,1.5,2];

% choose PAR range to look at
temp = -2:0.1:30;

%%%%%
%% nothing to be changed below this line
%%%%%

c10     = 10;
Tref    = 30;

% temp 
temp_lim = NaN(length(temp),length(Q10));
for i=1:length(Q10)
    temp_lim(:,i) = Q10(i).^(((temp + 273.15)-(Tref + 273.15))/c10);

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

for  i=1:numel(Q10)
    if ~exist('q10_string','var')
        q10_string = {strcat(['Q10=',num2str(Q10(i))])};
    else
        q10_string = [q10_string,strcat(['Q10=',num2str(Q10(i))])];
    end
end
clear i


fig = figure('Units','Centimeters', 'Position',[0 0 width height], 'Visible', vis);
for i=1:length(Q10)
    plot(temp,temp_lim(:,i),colors{i},'LineWidth',lw)
    hold on;
end
grid on;
axis([min(temp) max(temp) 0 1])
legend(q10_string,'Location','SouthEast')
xlabel('Temperature [°C]','FontSize',fs)
ylabel('temp limitation [n.d.]','FontSize',fs)
set(gca,'FontSize',fs)


