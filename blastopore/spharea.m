function SA = spharea(xx,yy,x0,y0,a)
% function SA = spharea(xx,yy,x0,y0,a)
%
% Calculate area on a sphere surface corresponding to projected
% area enclosed by a polygon in 2D, given 2D center and radius of
% the sphere and 2D polygon vertex coordinates. Uses Mapping Toolbox
% function areaint to calculate polygon projection surface area on
% sphere of radius a.
%
% ARGS:
% xx,yy = polygon vertices 
% x0,y0 = sphere center in (x,y)
% a     = sphere radius
%
% RETURNS:
% SA    = area on sphere surface corresponding to polygon projection
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 04/25/2003 JMT See JMT Notebook #2 pp117-118
%          07/15/2003 JMT Correct az/el transposition
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

% Operational flags
doplots = 0;

% Transform to sphere frame
xx = xx - x0;
yy = yy - y0;

% Calculate z coordinates of projection
a2 = a * a;
zz = sqrt(a2 - (xx.*xx + yy.*yy));

% Convert to spherical coordinates
[az,el,r] = cart2sph(xx,yy,zz);

% Convert to lat, lon
lat = el * 180 / pi;
lon = az * 180 / pi;

% Call Mapping Toolbox function which returns SA within
% polygon on S2 as a fraction of the total SA of the sphere.
% Post multiply by total SA of sphere of radius a
SA = areaint(lat,lon) * 4 * pi * a2;