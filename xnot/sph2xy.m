function [x,y] = sph2xy(la,lo,x0,y0,a0)
% [x,y] = sph2xy(la,lo,x0,y0,a0)
%
% Inverse transform for xy2sph
%
% ARGS :
% la,lo = latitude and longitude of (x,y) in the BFS frame
% x0,y0,a0 = BFC parameters
%
% RETURNS :
% x,y = vectors of coordinate components in 2D projection space
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/21/2003 JMT From scratch
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

la = la * pi / 180;
lo = lo * pi / 180;

[x,y,z] = sph2cart(lo,la,a0);

x = x + x0;
y = y + y0;