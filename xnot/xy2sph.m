function [la,lo] = xy2sph(x,y,x0,y0,a0)
% [la,lo] = xy2sph(x,y,x0,y0,a0)
%
% Convert an (x,y) pair from a 2D image of the sphere
% to lat and lon using the best fit sphere (circle) to
% define the spherical coordinate system.
%
% ARGS :
% x,y = vectors of coordinate components
% x0,y0,a0 = BFC parameters
%
% RETURNS :
% la,lo = latitude and longitude of (x,y) in the BFS frame
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

x = x - x0;
y = y - y0;

z = sqrt(a0.^2 - (x.^2 + y.^2));

[lo,la,r] = cart2sph(x,y,z);

lo = lo * 180 / pi;
la = la * 180 / pi;