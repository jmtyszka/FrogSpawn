function theta = dangle(O,A,B)
% theta = dangle(x0, y0, xa, ya, xb, yb)
%
% Calculate angular difference between vectors OA and OB
% by rotation to the frame which places OA on the +ve x axis.
%
% ARGS:
% O = origin column vector
% A = column vector A
% B = column vector B
%
% RETURNS:
% theta = angle AOB in degrees in range 0 to 360 degrees measured anticlockwise
%         from OA.
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/10/2003 JMT From scratch
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

% Vector differences
OA = A - O;
OB = B - O;

% Calculate 4-quadrant angles relative to x-axis using atan2
aOA = atan2(OA(2),OA(1)) * 180 / pi;
aOB = atan2(OB(2),OB(1)) * 180 / pi;

% Find difference and map to 0 to 360 domain
theta = mod(aOB - aOA,360);