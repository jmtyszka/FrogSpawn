function theta = dotangle(O,A,B)
% theta = dotangle(x0, y0, xa, ya, xb, yb)
%
% Calculate angular difference between vectors OA and OB
% from the scalar product relation
%    p = |a||b|cos(theta)
% => theta = acos(p / |a||b|)
%
% ARGS:
% O = origin column vector
% A = column vector A
% B = column vector B
%
% RETURNS:
% theta = angle AOB in degrees
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 06/19/2003 JMT From scratch
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

OA = A - O;
OB = B - O;

a = norm(OA);
b = norm(OB);

p = dot(OA,OB);

theta = acos(p / (a * b)) * 180 / pi; 