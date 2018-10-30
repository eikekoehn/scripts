%% author: Eike Koehn
% date Oct 17, 2018

%% set current path
cd /home/koehne/Documents/scripts/matlab

%% set MATLAB toolbox paths
addpath(genpath('/home/koehne/roms/roms_tools/'))

%% set data and grid paths
datapath = '/net/kryo/work2/fana/hindcast_bec/analysis/';
fdname = 'avg_O2_1979-2016_hindcast.nc';
datafile = [datapath,fdname];

gridpath = '/net/nardus/work/Ana/masks_and_maps/humpac_grids/';
fgname = 'humpac15_sm5rf015_grd.nc';
gridfile = [gridpath,fgname];

%% make a simple plot of an oxygen profile off the Peruvian coast
o2 = ncvarget(datafile,'O2',[100,100,1,1],[1,1,41,10]);

%%
figure
contourf(o2)
colorbar