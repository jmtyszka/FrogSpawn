function egg(eggdir)
% egg(eggdir)
%
% Wrapper function which calls current best fit egg shell function
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 02/05/2002 From scratch
%
% Example:
%   > egg;
% Runs best-fit sphere function in the current directory
%   > egg('data/embryos/20040613_SIM');
% Runs the function in 'data/embryos/20040613_SIM'
%
% See also eggprep, eggdmz, eggfan, eggshemi, eggshowbfs
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

if nargin < 1 eggdir = pwd; end

% Best fit sphere egg shell
egg_bfs(eggdir);