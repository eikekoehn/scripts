close all
clear all
clc

%% Add Data Path
path_name = '/net/kryo/work2/fana/hindcast_bec/hindcast/output/slavg/';
addpath(path_name)
disp('Data path added')

%% 
filename1 = 'slavg_1979.nc';
name = [path_name filename1];
ncid = netcdf.open(name,'NC_NOWRITE')

