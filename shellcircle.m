function [la, lo] = shellcircle(apole, maxlat, degs_samp)
% [la, lo] = shellcircle(apole, maxlat, degs_samp)
%
% Transform a uniform lat,lon sampling grid about an
% arbitrary pole back to the root angular coordinate system.
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 01/23/2002 From scratch
%
% Copyright 2002-2005 California Institute of Technology.
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

% Uniform sampling from the pole up to a maximum latitude
% in the pole frame of reference
la_0 = 0:degs_samp:maxlat;
lo_0 = 0:degs_samp:360;

[la_0_m,lo_0_m] = meshgrid(la_0,lo_0);

% Inverse transform to the root coordinate system
[la,lo] = rotatem(la_0_m,lo_0_m,apole,'inverse','degrees');