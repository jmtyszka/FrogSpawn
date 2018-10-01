function x = clamp(x,xlo,xhi)
% x = clamp(x,xlo,xhi)
%
% Restrict values of x to inclusive range [xlo,xhi]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 03/06/2002 From scratch (about time)
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

x(x>xhi) = xhi;
x(x<xlo) = xlo;