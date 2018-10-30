function  [zdata,sz] = zslice_matt(sdata,z,h,zdim,zeta, theta_s,theta_b,hc,N,type)
%ZSLICE: Interpolated ROMS nd-array sdata from sigma to z-levels
%  USAGE: 
%      [zdata,sz] = zslice(sdata,z,h,zdim,zeta,theta_s,theta_b,hc,N,type)
%    or
%      zdata = zslice(sdata,z,sz[,zdim])
%    with
%         sdata:  nd-array on sigma levels 
%             z:  target output depths vector.
%            sz: depth of sigma levels (output of zlevs)
%          zdim: sigma dimension (default: 3) i.e., sdata(x,y,z,....)
%           h, zeta, theta_s, theta_b, hc, N, type: as in zlevs.m
%       
%  CRUDE CHECK THAT DIMENSION OK
size1h = size(h,1);
size1d = size(sdata,1);
size2h = size(h,2);
size2d = size(sdata,2);
if length(size(sdata))==3
if size1h~=size1d || size2h~=size2d;  display('Check dimensions!');stop;end
elseif length(size(sdata))==2
    if size1h~=size1d ;  display('Check dimensions!');stop;end
else
    stop
end

if nargin<4
    zdim = 3;
end
if zdim > 1
    % pidx: index to permute sigma-coord to 1st dimension:
    nd = length(size(sdata));
    zidx=find(size(sdata)==N);
    pidx = [zidx,1:zidx-1,zidx+1:nd];
    %pidx = [zdim,1:zdim-1,zdim+1:nd];
    sdata = permute(sdata,pidx);
end
if nargin>5
    if nargin<10
        error('roms:zslice',...
            'Not enough input variables to compute depth of sigma levels');
    end
    if nargin < 5
        zeta = 0;
    end
    sz = squeeze(zlevs(h,zeta,theta_s,theta_b,hc,N,type));
end
mex /net/bellis/work/ifrenger/ROMS_input/Roms_tools/Preprocessing_tools/interp1nd.c ;
system(['mv interp1nd.mexa64 ./']);
zdata = interp1nd(sz,squeeze(sdata),z');
if zdim > 1 
    % permute back
    zdata = ipermute(zdata,pidx);
end
return
end
