function tissottest
% tissottest
%
% Draw tissot circles in orthographic projection of the south pole to
% demonstrate spherical distortion.
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/20/2004 JMT From scratch
%
% Copyright 2004-2005 California Institute of Technology.
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

figure(1); clf; colormap(gray)

load topo

topo = ones(size(topo));

axesm('mapprojection','ortho','origin',[-90 0 0]);

meshm(topo,topolegend);

tissot([30 60 0.2],'color','w','linewidth',2);

% Add circle at south pole
[la,lo] = scircle1(-89.9,0,0.2*180/pi)
plotm(la,lo,'w','linewidth',2);
