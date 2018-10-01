function res = xnotanal(xx,yy,x0,y0,a0,x_AB,y_AB)
% res = xnotanal(xx,yy,x0,y0,a0,x_AB,y_AB)
%
% Geometric analysis of the XNot domain
% - projected and surface area ratios to total embryo
% - min, max, mean, median width of XNot domain
% 
% ARGS:
% xx,yy    = XNot domain boundary vertices (flat projection)
% x0,y0,a0 = best fit circle parameters (center and radius)
% x_A,y_A  = XNot domain long axis GC definition (two points)
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/16/2003 JMT From scratch
%
% Copyright 2003-2005 California Institute of Technology.
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

% Define BFS geoid
bfs_geoid = [a0 0];

% Define angular sampling resolution in degrees
asamp = 0.5;

%-----------------------------------------------------------------------
% Analyze the XNot domain area
%-----------------------------------------------------------------------

% Calculate projected and S2 areas and ratios
res.xnot_area     = polyarea(xx,yy);
res.xnot_surfarea = spharea(xx,yy,x0,y0,a0);
res.emb_area      = pi * a0 * a0;
res.emb_surfarea  = 2 * res.emb_area;

%-----------------------------------------------------------------------
% Rotate the XNot domain to the equator using the axis GC as a guide
%-----------------------------------------------------------------------

% Transform all coordinates to the BFS spherical frame
[la_x,lo_x] = xy2sph(xx,yy,x0,y0,a0);
[la_AB,lo_AB] = xy2sph(x_AB,y_AB,x0,y0,a0);

% Extract axis GC components
la_A = la_AB(1); lo_A = lo_AB(1);
la_B = la_AB(2); lo_B = lo_AB(2);

%-----------------------------------------------------------------------
% Rotate coordinate system
%-----------------------------------------------------------------------

% Find normal to plane through AOB
% See JMT Notebook #2 pp156

% Normal vector Euler rotations
[Ax,Ay,Az] = sph2cart(lo_A*pi/180,la_A*pi/180,1);
[Bx,By,Bz] = sph2cart(lo_B*pi/180,la_B*pi/180,1);

AOBn = cross([Ax Ay Az],[Bx By Bz]);

fprintf('Normal to AOB : (%0.3f, %0.3f, %0.3f)\n', AOBn(1), AOBn(2), AOBn(3));

[THn,PHIn,Rn] = cart2sph(AOBn(1),AOBn(2),AOBn(3));

fprintf('Normal to AOB : Az %0.3f deg  El %0.3f deg\n', THn * 180/pi, PHIn * 180/pi);

alpha = -THn;
beta  = -pi/2-PHIn;

% Rotation angles
ca = cos(alpha); sa = sin(alpha);
cb = cos(beta); sb = sin(beta);

% Construct Euler rotation matrices
Rz = [ca -sa 0; sa ca 0; 0 0 1];
Ry = [cb 0 -sb; 0 1 0; sb 0 cb];

% Convert XNot to cartesian
[xx,xy,xz] = sph2cart(lo_x*pi/180,la_x*pi/180,1);

% Flatten into 3 x n matrix
xnot_c = [xx(:)';xy(:)';xz(:)'];

% Euler rotations
xnot_eq = Ry * Rz * xnot_c;

% Extract coordinates
xx_eq = xnot_eq(1,:);
xy_eq = xnot_eq(2,:);
xz_eq = xnot_eq(3,:);

% Convert back to spherical coordinates
[lo_x0,la_x0,r_x0] = cart2sph(xx_eq,xy_eq,xz_eq);

la_x0 = la_x0 * 180/pi;
lo_x0 = lo_x0 * 180/pi;

%-----------------------------------------------------------------------
% Size analysis of XNot domain on equator
%-----------------------------------------------------------------------

% Find XNot boundary box to nearest degree
min_lax = floor(min(la_x0(:)))-1;
max_lax = ceil(max(la_x0(:)))+1;
min_lox = floor(min(lo_x0(:)))-1;
max_lox = ceil(max(lo_x0(:)))+1;

% Create a coordinate mesh within the boundary box
la_v = min_lax:asamp:max_lax;
lo_v = min_lox:asamp:max_lox;
[lo_m,la_m] = ndgrid(lo_v,la_v);

% lo_m and la_m are nlo x nla in size

% Find interior of XNot domain
fprintf('Finding XNot domain interior points\n');
[inx,onx] = inpolygon(lo_m,la_m,lo_x0,la_x0);

% Add boundary points
inx = inx | onx;

figure(2); clf; colormap(gray);
hold on;
plot(lo_x0,la_x0);
plot(lo_m(inx),la_m(inx),'r.');
axis equal
set(gca,'XDir','reverse');
xlabel('Longitude');
ylabel('Latitude');
hold off;

% Loop over longitudes recording latitude limit stats

fprintf('Analyzing XNot domain angular widths\n');

% Note dimension order from ndgrid above
[nlo,nla] = size(lo_m);
la_width = zeros(nlo,1);

for lc = 1:nlo

  ii = find(inx(lc,:));
  if ~isempty(ii)
    la_width(lc) = la_v(max(ii)) - la_v(min(ii));
  else
    la_width(lc) = 0;
  end
  
end

% Remove zeros from la_width
la_width(la_width == 0) = [];

res.x_len = max_lox - min_lox;
res.y_min = min(la_width(:));
res.y_max = max(la_width(:));
res.y_mean = mean(la_width(:));
res.y_med = median(la_width(:));

res.x_y = res.x_len / res.y_med;
