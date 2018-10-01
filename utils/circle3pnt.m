function [r,x0,y0] = circle3pnt(x,y)
% [r,x0,y0] = circle3pnt(x,y)
%
% Calculate circle center and radius in 2D from 
% three 2D points.
%
% ARGS:
% x = three x coords (3 x 1)
% y = three y coords (3 x 1)
%
% RETURNS:
% r  = circle radius
% x0 = x coord of circle center
% y0 = y coord of circle center
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 4/21/2003 JMT From scratch - Use Mathematica solution
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

% Default returns
r = 0.0;
x0 = 0.0;
y0 = 0.0;

if length(x) < 3 return; end
if length(y) < 3 return; end

% Extract elements
x1 = x(1); y1 = y(1);
x2 = x(2); y2 = y(2);
x3 = x(3); y3 = y(3);

% Calculate x0 (See Mathematica soln)
x0 = ((x3*x3) * (y1 - y2) + ...
  (x1*x1 + (y1 - y2)*(y1 - y3))*(y2 - y3) + (x2*x2) *(-y1 + y3)) / ...
  (2*(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3)));
  
y0 = -(-2*(x2 - x3) * (-x1*x1 + x3*x3 - y1*y1 + y3*y3) + ...
  2*(x1 - x3)*(-x2*x2 + x3*x3 - y2*y2 + y3*y3)) / ...
  (4.*(x3*(y1 - y2) + x1*(y2 - y3) + x2*(-y1 + y3)));

% Substitute into equation for r
r = sqrt((x1-x0)^2 + (y1-y0)^2);     