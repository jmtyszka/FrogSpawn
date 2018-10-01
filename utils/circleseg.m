function [xx,yy] = circleseg(x0,y0,a,segrng,varargin)
% [xx,yy] = circle(x0,y0,a,segrng,plotopts...)
%
% Draw a circle in an axis
%
% ARGS:
% x0,y0    = circle center
% a        = circle radius
% segrng   = segment range in degrees (eg [-90 30])
% varargin = standard plot options (eg 'r:','linewidth',3,...)
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 10/19/2001 JMT from scratch
%          05/16/2003 JMT Add return arguments
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

% Default args
if nargin < 1 x0 = 0; end
if nargin < 2 y0 = 0; end
if nargin < 3 R = 1; end
if nargin < 4 segrng = [-180 180]; end

% Define circle at 1 degree intervals
dtheta = 1;
theta = (segrng(1):dtheta:segrng(2)) * pi / 180;

% Vertices
xx = a * cos(theta) + x0;
yy = a * sin(theta) + y0;

% Plot circle if no return arguments requested
if nargout < 1
  plot(xx,yy,varargin{:});
end