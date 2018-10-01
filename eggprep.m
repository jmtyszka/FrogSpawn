function eggprep(datadir,ndig,vsize)
% eggprep(datadir,ndig,vsize)
%
% Prepare a high-resolution egg dataset for processing.
% Assume a name format of <name>xxxx.tif
%
% ARGS:
% datadir = data directory containing tif files
% ndig    = number of zero padded digits in image number [4]
% vsize   = nominal voxel size (dx dy dz) [1 1 1]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 08/02/2001 JMT From scratch
%          08/22/2001 JMT Generalize for use with other egg data
%          05/02/2003 JMT Add ndigit argument
%          03/16/2004 JMT Remove calculation of voxel center grid due to
%                         maximum variable size limits
%          09/22/2004 JMT Make interactive directory choice optional
%          11/02/2004 JMT Add arguments allowing batch operation
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

% Software version
fver = 3.3;

% Default directory
datadir = pwd;

% Prompt for an example filename from the series
[fname, dname] = uigetfile(fullfile(datadir,'*.tif'), 'Select one TIF from the series');

% Default arguments
if nargin < 2 ndig = 3; end
if nargin < 3 vsize = [1 1 1]; end

fprintf('Using a voxel size of (%0.1f x %0.1f x %0.1f)um\n', vsize(1), vsize(2), vsize(3));

% Catch cancel button
if isequal(fname,0) | isequal(dname,0)
  fprintf('No file selected - exiting\n');
  return
end

% Get information from the selected TIF image
info = imfinfo(fullfile(dname,fname));

% Check for TIF format image file
if ~strcmp(info.Format,'tif')
  fprintf('Selected image is not a TIF file - exiting\n');
  return
end

nx = info.Height;
ny = info.Width;

% Extract the stub (all but last eight characters)
fstub = fname(1:(length(fname)-4-ndig));

% Create a format string for the image filename
fstr = sprintf('%%0%dd.tif',ndig);

% Loop upwards through image indices to find maximum index
fprintf('Searching for TIF images with this name format\n');
nz = 1;
imname = fullfile(dname,[fstub sprintf(fstr,nz)]);
while exist(imname)
  nz = nz + 1;
  imname = fullfile(dname,[fstub sprintf(fstr,nz)]);
end

% Bail out if no images found
if nz < 1
  fprintf('No images with format %s found\n', fstr);
  return
end
  
nz = nz - 1;
fprintf('  %d images located\n', nz);

% Allocate uint8 memory
fprintf('Allocating %d MB (%d x %d x %d)\n', round(nx * ny * nz / 1024^2), nx, ny, nz);
s = zeros8([nx ny nz]);

% Load data slice by slice as uint8
fprintf('Loading images ');
zstep = fix(nz/10);
for z = 1:nz
  if mod(z,zstep) == 0
    fprintf('-');
  end
  imname = fullfile(dname,[fstub sprintf(fstr,z)]);
  this_s = uint8(imread(imname));
  s(:,:,z) = squeeze(this_s(:,:,1));
end
fprintf(' done\n');

% Save data
matname = fullfile(dname,[fstub '.mat']);
fprintf('Saving egg image data in %s\n', matname);
save(matname,'s','vsize');