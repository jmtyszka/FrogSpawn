function [r_fit,CofM] = fitsphere(x,y,z)
% [r_fit,CofM] = fitsphere(x,y,z)
%
% Non-linear least-squares fit of an sphere to (x,y,z) triples
%
% ARGS:
% x,y,z = coordinate vectors of points on sphere (N x 1)
%
% RETURNS:
% CofM  = center of sphere vector [x0,y0,z0] in same units as x,y,z
% r_fit = sphere radius in same units as x,y,z
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/30/2001 See Eberly, D Magic Software article
%          01/16/2002 Change output argument style for egg.m
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

% Operation flags
showbfs = 0;

% Constants
max_count = 1000;

% Flatten the coordinate matrices
x = x(:);
y = y(:);
z = z(:);

% Calculate mean locations of points
x_m = mean(x);
y_m = mean(y);
z_m = mean(z);

% Initial guess
a = x_m;
b = y_m;
c = z_m;

converged = 0;
count = 0;

while ~converged
  
  % Calculate L_i, etc
  L_i = sqrt((x-a).^2 + (y-b).^2 + (z-c).^2);
  L_m = mean(L_i); 
  L_a = mean((a-x)./L_i);
  L_b = mean((b-y)./L_i);
  L_c = mean((c-z)./L_i);

  % Calculate new estimate of sphere center
  a_n = x_m + L_m * L_a;
  b_n = y_m + L_m * L_b;
  c_n = z_m + L_m * L_c;
  
  % Calculate difference
  d = norm([(a_n - a) (b_n - b) (c_n - c)]) / norm([a b c]);
  
  % Increment counter
  count = count + 1;

  if (d < 1e-3) | (count > max_count)
    converged = 1;
  else
    a = a_n;
    b = b_n;
    c = c_n;
  end
  
end

% Pass return arguments
CofM  = [a b c];
r_fit = L_m;

if showbfs
  
  figure(6); clf
  set(gcf,'color','k');
  
  % Plot the sample vertices
  plot3(x,y,z,'g'); hold on;

  % Scale unit sphere to BFS
  [sx,sy,sz] = sphere(128);
  sx = sx * r_fit + a;
  sy = sy * r_fit + b;
  sz = sz * r_fit + c;
  
  % Plot the BFS
  surf(sx,sy,sz); hold off;
  
end



