function blast(imstep)
% blast(imstep)
%
% Use closed natural spline (cscvn) to specify blastopore and embryo outlines
% on an image under user guidance.
% - UI allows addition and deletion from points in the spline.
% - Enclosed areas calculated from points list generated from spline
%   function (fnplt).
% - Spherical correction applied to blastopore area.
% - Expects tif stacks (see ImageJ) as input
%
% ARGS:
% imstep = image step size through stack [5]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 04/19/2003 JMT From scratch
%          05/15/2003 JMT Add support for stacks
%          07/10/2003 JMT Tidy up args and output
%          07/15/2003 JMT Correct SA calculation in blastproj.m (JMT)
%          08/07/2003 JMT Change default imstep to 3 (as per AJE)
%          08/07/2003 JMT Add output only for BP closure (1 - area ratio)
%          03/15/2004 JMT Thicken plot lines, enlarge plot markers and
%                         replot BFS points
%          10/16/2004 JMT Add support for color eggpop rendering
%                         Correct axis orientation in eggpop
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

% Default parameters
if nargin < 1 imstep = 3; end

% Ask user for example image from series
[fname,dname] = uigetfile('*.tif');

if isequal(fname,0) | isequal(dname,0)
  return
end

% Construct full path to image file
imname = fullfile(dname,fname);

% Get image info
try
  fprintf('Loading TIFF stack\n');
  info = imfinfo(imname);
catch
  fprintf('Unexpected error reading image information\n');
  return;
end

% Parse image info
nims = length(info);

%----------------------------------------------------
% Open a results file
%----------------------------------------------------
fstub = fname(1:(length(fname)-4));
resname = fullfile(dname,[fstub '_res.txt']);
fd = fopen(resname,'w');
if fd < 1
  fprintf('Problem opening results file\n');
  return
end

% Column headers for results file
fprintf(fd, 'Im_No Area_Ratio Blast_SA Emb_SA Blast_A Emb_A x0 y0 a0 x_pnts y_pnts\n');

% Initialize BFS parameters
x0 = [];
y0 = [];
a0 = [];

% Create a new grayscale figure
figure(1); clf; colormap(gray);
set(gcf,'NumberTitle','off');
set(gcf,'Units','Centimeters','Position',[1 1 15 15]);
set(gcf,'Renderer','zbuffer');

%----------------------------------------------------
% Image loop over all frames of movie
%----------------------------------------------------

for ic = 1:imstep:nims 
  
  % Load image
  s = imread(fname,ic);
  
  % Draw initial image
  figure(1); imagesc(s); axis equal off; title(sprintf('%s : %d/%d',fname,ic,nims));
  
  %----------------------------------------------------
  % Define embryo border
  %----------------------------------------------------
  
  % Initialize reuse BFS flag
  reuse_bfs = 0;
  
  % Display previous BFS
  if ~isempty(a0)
    
    % Overlay the circle
    hold on;
    circle(x0,y0,a0,'g','linewidth',3);
    hold off;
    
    % Prompt user in figure title bar
    set(gcf,'Name','Press SPACE to accept previous BFS. Q to quit');
    
    % SPACE accepts, anything else rejects
    if waitforbuttonpress == 1
      switch get(gcf,'CurrentCharacter')
        case ' '
          reuse_bfs = 1;
        case 'q'
          fprintf('Quitting at user request\n');
          fclose(fd);
          return
      end
    end
    
  end
  
  % Redefine BFS if necessary
  if reuse_bfs == 0
    
    keep_going = 1;
    
    while keep_going
      
      % Redraw initial image to clear any previous BFS
      figure(1);
      imagesc(s); axis equal off; title(sprintf('%s : %d/%d',fname,ic,nims));
      
      set(gcf,'Name','Select three points on embryo boundary');
      
      % Initialize point vectors
      xbfc = zeros(3,1);
      ybfc = zeros(3,1);
      
      for p = 1:3
        
        % Get (x,y) coordinate from user
        [xbfc(p),ybfc(p),btn] = ginput(1);
        
        hold on;
        plot(xbfc,ybfc,'og','markersize',12,'linewidth',3);
        hold off;
        
      end
      
      % Calculate specified circle from the three points
      [a0,x0,y0] = circle3pnt(xbfc,ybfc);
      
      % Overlay the circle
      hold on;
      circle(x0,y0,a0,'g','linewidth',3);
      hold off;
      
      % Ask user to confirm BFC
      ans = questdlg('Confirmation','Accept Best Fit Circle?','Yes','No','Yes');
      
      switch ans
        case 'Yes'
          keep_going = 0;
        case 'No'
          keep_going = 1;
      end
      
    end
    
  end
  
  %----------------------------------------------------
  % Define blastopore border
  %----------------------------------------------------
  
  set(gcf,'Name','Define blastopore boundary. Backspace deletes point. Space accepts all. Q quits');
  
  % Initialize continuation flag
  keep_going = 1;
  
  % Initialize cumulative points lists and counter
  n = 0;
  pnts = [0;0];
  
  while keep_going
    
    %----------------------------------------------------
    % Get (x,y) coordinate from user
    %----------------------------------------------------
    [x,y,button] = ginput(1);
    
    %----------------------------------------------------
    % Handle button and key presses
    %----------------------------------------------------
    switch button
      
      case 8 % Delete key
        
        % Eliminate previous point
        if n > 0
          pnts(:,n) = [];
          n = n - 1;
        end
        
      case 32 % Space key completes
        
        keep_going = 0;
        
      case {81,113} % q,Q key
        
        fprintf('Quitting at user request\n');
        fclose(fd);
        return
        
      otherwise
        
        % Overwrite last point with new point
        n = n + 1;
        pnts(:,n) = [x;y];
        
    end
    
    if keep_going
      
      % Refresh the figure
      figure(1);
      imagesc(s); axis equal off; title(sprintf('%s %d/%d',fname,ic,nims));
      
      % Redraw the BFC
      hold on;
      plot(xbfc,ybfc,'og','markersize',12,'linewidth',3);
      circle(x0,y0,a0,'g','linewidth',3);
      hold off;
      
      %----------------------------------------------------
      % Calculate closed natural spline through all points
      %----------------------------------------------------
      if n > 0
        
        cs = cscvn([pnts pnts(:,1)]);
        
        % Calculate spline plot points
        pp = fnplt(cs);
        xx = pp(1,:); yy = pp(2,:);
        
        %----------------------------------------------------
        % Overplot spline on image
        %----------------------------------------------------
        
        % Redraw the spline
        hold on
        plot(xx,yy,'r','linewidth',3);
        plot(pnts(1,:),pnts(2,:),'or','markersize',12,'linewidth',3);
        hold off
        
      end
      
    end
    
  end
  
  if isempty(xx) | isempty(yy)
    
    bname = questdlg('Blastopore boundary unspecified',...
      'Repeat blastopore boundary?',...
      'Repeat','Continue');
    
  end
  
  % Calculate blastopore surface area from projected area enclosed by
  % spline. See JMT Notebook #2 p118
  fprintf('Calculating area ratios from %d line samples\n', length(xx));
  blast_area = polyarea(xx,yy);
  blast_surfarea = spharea(xx,yy,x0,y0,a0);
  emb_area = pi * a0 * a0;
  emb_surfarea = 2 * emb_area;
  
  proj_area_ratio = blast_area / emb_area;
  proj_surfarea_ratio = blast_surfarea / emb_surfarea;
  
  % Calculate closures (1 - area)
  BP_closure_area = 1 - proj_area_ratio;
  BP_closure_surfarea = 1 - proj_surfarea_ratio;
  
  fprintf('------------------------------------------\n');
  fprintf('Filename : %s\n',fname);
  fprintf('Image    : %d\n', ic);
  fprintf('------------------------------------------\n');
  fprintf('Blastopore:\n');
  fprintf('  Projected area     : %0.5f\n', blast_area);
  fprintf('  Surface area       : %0.5f\n', blast_surfarea);
  fprintf('Embryo:\n');
  fprintf('  Projected area     : %0.5f\n', emb_area);
  fprintf('  Surface area       : %0.5f\n', emb_surfarea);
  fprintf('------------------------------------------\n');
  fprintf('Projected Area Ratio : %0.5f\n', blast_area / emb_area);
  fprintf('Surface Area Ratio   : %0.5f\n', blast_surfarea / emb_surfarea);
  fprintf('------------------------------------------\n');
  fprintf('BP Closure by Area   : %0.5f\n', BP_closure_area);
  fprintf('BP Closure by SA     : %0.5f\n', BP_closure_surfarea);
  fprintf('------------------------------------------\n');
  fprintf('\n');
  
  % Append results to file
  fprintf(fd, '%d ', ic);
  fprintf(fd, '%0.5f ', blast_surfarea / emb_surfarea);
  fprintf(fd, '%0.5f ', blast_surfarea);
  fprintf(fd, '%0.5f ', emb_surfarea);
  fprintf(fd, '%0.5f ', blast_area);
  fprintf(fd, '%0.5f ', emb_area);
  fprintf(fd, '%0.5f %0.5f %0.5f ', x0, y0, a0);
  fprintf(fd, '%0.5f ', pnts(1,:)); % x points for blastopore 
  fprintf(fd, '%0.5f ', pnts(2,:)); % y points for blastopore
  fprintf(fd, '\n');

end % Image loop

% Close results file
fclose(fd);

%--------------------------------------------------------------
% Create a surface mesh for the image
%--------------------------------------------------------------

% Note axis swap - x is the second coordinate (columns)
[ny,nx,nc] = size(s);
[xc,yc] = meshgrid((1:nx)-x0-1,(1:ny)-y0-1);
r2 = xc.*xc + yc.*yc;
a2 = a0 * a0;
r2(r2 > a2) = a2;
zc = sqrt(a2 - r2);

%--------------------------------------------------------------
% Render the embryo surface in 3D
%--------------------------------------------------------------

figure(2); clf;
h = surface(xc,yc,zc);
set(h,'CData',s,'FaceColor','texturemap','EdgeColor','none');

axis equal vis3d ij off;
camproj('perspective');

set(gcf,'color',[0.8 0.8 0.8]);

% Draw blastopore boundary
xx = xx - x0;
yy = yy - y0;
zz = sqrt(a2 - (xx.^2 + yy.^2));
hold on
plot3(xx,yy,zz,'r-','linewidth',3);
hold off