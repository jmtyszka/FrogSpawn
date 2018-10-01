function egg_bfs(eggdir)
% egg_bfs(eggdir)
%
% Extract spherical shells for surface mapping of frog embryos using the best fit
% sphere (BFS) approximation.
%
% ARGS:
% eggdir = optional egg dataset directory
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/30/2001 0.1 From scratch
%          08/02/2001 0.2 Add support for full resolution data with downsampling
%                         for sphere fit
%                         Add mapping support
%          08/15/2001 0.3 Add test shell for south pole specification
%                         Add multiple rscale values with loop to reduce out-of-memory errors
%          08/16/2001 0.4 Add south pole transform
%          08/21/2001 0.5 Correct polar rotation transform
%          01/22/2002 0.6 Complete upgrade to spherical harmonics
%          01/30/2002 1.0 Split development into BFS and SH fits
%          02/05/2002 1.0 Remove unused arguments and returns
%          03/06/2002 1.1 Replace blastopore test shell with 3-plane interactive view
%          05/30/2002 1.2 Improve robustness of binary segmentation
%          -------------- Skip version 2.0 for consistency with FrogSpawn
%                         versioning
%          03/16/2004 3.1 Calculate voxel center grid on the fly
%
% Copyright 2001-2005 California Institute of Technology.
% All rights reserved.
%
% This file is part of FrogSpawn.
%
% FrogSpawn is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% FrogSpawn is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with FrogSpawn; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% Various parameters
eggver = 3.1;                  % Software version number
meshfactor = 0.1;              % Mesh reduction factor for test isosurface
rscale_bptest = [0.7 0.8 0.9]; % Fractional radii for blastopore test shells
ds_bptest = 2;                 % Degrees per sample for the blastopore test shell
ds_full = 0.15;                % Degrees per sample for the full res shell
allpnts = [-90 90 -180 180];   % Lat and Lon limits for whole shell

% Operation flags
findae     = 1;         % Locate archenteron
hardcpy    = 1;         % Hardcopy important figures to files
savemesh   = 0;         % Save the isosurface mesh
saveinside = 1;         % Save the "inside egg" mask

% Defaults
if nargin < 1
  eggdir = 'C:\Shared Data\Frog Eggs\Frog739';
end

% Splash text
fprintf('\n');
fprintf('EGG CARTOGRAPHY - BEST FIT SPHERE\n');
fprintf('---------------------------------\n');
fprintf('Version     : %0.1f\n', eggver);
fprintf('Copyright   : 2002 California Institute of Technology\n');
fprintf('Programming : Mike Tyszka\n');
fprintf('Design      : Mike Tyszka and Andy Ewald\n');
fprintf('\n');

% Allow user to select the egg data file from datadir
[eggname, eggdir] = uigetfile(fullfile(eggdir,'*.mat'));
if isequal(eggname,0) | isequal(eggdir,0)
  fprintf('No egg data selected - exiting\n');
  return
end

% Extract file stub
eggstub = eggname(1:(length(eggname)-4));

% Load uint8 egg data (s) and voxel dimension vector (vsize)
fprintf('Loading data from %s\n', eggname);
s = []; % Default
eggmat = fullfile(eggdir,eggname);
load(eggmat);

% Check for existance of s[]
if isempty(s)
  fprintf('No egg data (s matrix) in %s - exiting\n', eggname);
  return
end

% Grab dimensions
[nx,ny,nz] = size(s);

%----------------------------------------------------
% Conditional downsampling depends on voxel count
%----------------------------------------------------

% Create a downsampled version if necessary for sphere fitting
if prod([nx ny ny]) > 128^3
  dwnfac = 6;
else
  dwnfac = 1;
end

fprintf('  Downsampling by %d for segmentation\n', dwnfac);
s_dwn = reducevolume(s, [dwnfac dwnfac dwnfac]);
s_dwn = smooth3(s_dwn,'gaussian');
[nxd,nyd,nzd] = size(s_dwn);

%----------------------------------------------------
% Binary segment the egg
%----------------------------------------------------

fprintf('  Segmenting outside of egg\n');
s_in = seggment(s_dwn);

% Preview the surface segmentation
figure(1); clf; colormap(gray);
set(gcf,'Name','Segmentation preview');
s_m = bic_montage(s_dwn);
s_in_m = bic_montage(s_in);
s_m = s_m / max(s_m(:));
s_in_m = s_in_m / max(s_in_m(:));
[nxm,nym] = size(s_m);
s_rgb = zeros(nxm,nym,3);
s_rgb(:,:,1) = s_in_m;
s_rgb(:,:,2) = s_m;
imagesc(s_rgb);
axis equal off;
drawnow

%----------------------------------------------------
% Best fit sphere to downsampled data
%----------------------------------------------------

% Generate voxel center grid (um)
vsize_dwn = vsize * dwnfac;
xv = ((1:nxd)-0.5) * vsize_dwn(1);
yv = ((1:nyd)-0.5) * vsize_dwn(2);
zv = ((1:nzd)-0.5) * vsize_dwn(3);
[xcd,ycd,zcd] = ndgrid(xv,yv,zv);

% Generate an isosurface for the egg
% Voxel coordinate meshes based on original sampling grid (um)
fprintf('  Generating isosurface\n');
fv = isosurface(xcd, ycd, zcd, s_in, 0.5);

% Optional mesh save
if savemesh
  save(fullfile(eggdir,'eggmesh.mat'),'fv');
end

% Optional segmented data save
if saveinside
  save(fullfile(eggdir,'egginside.mat'),'s_dwn','s_in');
end

% Extract the vertices (um)
fprintf('  Extracting vertices : ');
x_samp = fv.vertices(:,1);
y_samp = fv.vertices(:,2);
z_samp = fv.vertices(:,3);
fprintf('%d\n', length(x_samp));

% Get the best fit sphere parameters
fprintf('  Fitting sphere to vertices\n');
[r0,CofM] = fitsphere(x_samp,y_samp,z_samp);

% Transform coordinate mesh to CofM frame (um)
xcd = xcd - CofM(1);
ycd = ycd - CofM(2);
zcd = zcd - CofM(3);

% Report BFS results
fprintf('BEST FIT SPHERE\n');
fprintf('  Radius : %0.1f um\n', r0);
fprintf('  Center : (%0.1f, %0.1f, %0.1f) um\n', CofM(1), CofM(2), CofM(3));

%---------------------------------------------------------------------------
% Draw three-plane view through center of mass with sphere section overlaid
%---------------------------------------------------------------------------

figure(2); clf; colormap(gray);
% set(gcf,'Name','Best Fit Sphere','Color','k','InvertHardCopy','off');

% Identify blastopore center in a 3-plane view
[xbp,ybp,zbp] = egg3plane(xcd, ycd, zcd, s_dwn, r0, [0 0 0]);

% Optional hardcopy of 3-plane BFS views in current figure
if hardcpy
  figure(2);
  print('-dtiff','-r300',fullfile(eggdir,[eggstub 'BFS.tif']));
end

% Convert blastopore center coordinates to lat and lon
[azbp,elbp,rbp] = cart2sph(xbp,ybp,zbp);
bp(1) = elbp * 180/pi;
bp(2) = azbp * 180/pi;

fprintf('Blastopore            : %0.1f deg Lat, %0.1f deg Long\n', bp(1), bp(2));

% Create a new mapping origin which places the blastopore at the south pole
spole = bp;
npole = antipode(spole(1),spole(2));
bporigin = newpole(npole(1),npole(2));

%----------------------------------------------------------------
% Locate the archenteron within several latitude planes
%----------------------------------------------------------------

if findae

  figure(4); clf; colormap(gray);
  set(gcf,'Name','Locate center of archenteron','color','k');
  
  % Latitudes to extract
  la_ae = -60:5:-20;

  % Extract latitude planes for archenteron location
  [lat_s, lat_x, lat_y, lat_z] = eggplane_bfs(uint8(s_dwn),vsize_dwn,la_ae,r0,CofM,npole,0,'lat');
  
  % Draw the latitude planes
  for lac = 1:9
    
    s_l = lat_s(:,:,lac);
    x_l = lat_x(:,:,lac);
    y_l = lat_y(:,:,lac);
    
    subplot(3,3,lac), pcolor(x_l, y_l, s_l);
    shading interp;
    axis equal xy off;
    title(sprintf('Latitude %0.1f',la_ae(lac)),'color','w');
    
  end

  % Let user locate the archenteron center
  ae_xy = ginput(1);
  
  % Translate to map coordinate system
  ae = atan2(ae_xy(2), ae_xy(1)) * 180 / pi;

  fprintf('Archenteron longitude : %0.1f deg\n', ae);
  
  % Optional hardcopy of AE finder in current figure
  if hardcpy
    figure(4);
    print('-dtiff','-r300',fullfile(eggdir,[eggstub 'AE.tif']));
  end

end

return

%-------------------------------------------------------
% Save best fit sphere, pole and archenteron parameters
%-------------------------------------------------------
eggstub = eggname(1:(length(eggname)-4));
parfile = [eggdir eggstub '_bfs.mat'];
fprintf('Saving BFS parameters in %s\n', parfile);
save(parfile,'r0','CofM','npole','ae', 'vsize');