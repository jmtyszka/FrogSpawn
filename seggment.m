function s_in = seggment(s)
% s_in = seggment(s)
%
% ARGS:
% s = 3D scalar field
%
% RETURNS:
% s_in = inside-outside function for egg
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 01/22/2002 JMT From scratch
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

% Operation flags
saveseg = 1;         % Save segmentation images

% Grab data dimensions
[nx,ny,nz] = size(s);

% Make space for segmented image
s_in = s;

% Use Otsu's Method to determine optimal grayscale threshold
th = graythresh(s);
bw = double(s >= th);
  
% Create a waitbar
hwb = waitbar(0,'');

% Loop over all XY slices
for z = 1:nz
  
  % Update waitbar
  waitbar(z/nz,hwb,'Segmenting XY slices');
  
  % Extract current XY slice
  bwz = bw(:,:,z);
  
  % Fill holes
  bw_noholes = imfill(bwz,'holes');
  
  % Morphological opening (removes small islands)
  se = strel('disk',3);
  bw_open = imopen(bw_noholes,se);
  
  % Store segmented image
  s_in(:,:,z) = double(bw_open);
  
  if saveseg & (z == fix(nz/2))
    sz = s(:,:,z);
    save eggsegment sz bw bw_noholes bw_open;
  end
  
end

% Close the waitbar
close(hwb);
