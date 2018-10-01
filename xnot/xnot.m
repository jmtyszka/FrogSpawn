function xnot(imstep)
% xnot(imstep)
%
% Use closed natural spline (cscvn) to specify XNot domain and embryo outlines
% on an image under user guidance.
% - UI allows addition and deletion from points in the spline.
% - Enclosed areas calculated from points list generated from spline
%   function (fnplt).
% - Spherical correction applied to Xnot area and dimensions.
% - Expects tif stacks (see ImageJ) as input
%
% ARGS:
% imstep = image step size through stack [5]
%
% AUTHOR : Mike Tyszka, Ph.D.
% PLACE  : Caltech BIC
% DATES  : 07/16/2003 JMT Clone from blast.m (JMT)
%          08/07/2003 JMT Fixed bug in matrix indexing due to ndgrid use
%          08/07/2003 AJE Commented out screen reporting of Area details
%          Proj Area Ratio and Min Width
%          03/15/2004 JMT Make overlay lines thicker and replot BFC points
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
if nargin < 1 imstep = 5; end

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
fprintf(fd, 'Im_No SA_Ratio L_MEDW XNot_SA Emb_SA Xnot_A Emb_A Xlen Ymin Ymax Ymean Ymed x0 y0 a0 XNot_x_pnts XNot_y_pnts XNot_axis_lats XNot_axis_lons\n');

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
  % Define XNot domain boundary
  %----------------------------------------------------
  
  set(gcf,'Name','Define XNot domain. Backspace deletes point. Space accepts all. Q quits');
  
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
    
    bname = questdlg('XNot domain unspecified',...
      'Repeat XNot domain?',...
      'Repeat','Continue');
    
  end
  
  %----------------------------------------------------
  % User defined long axis of XNot domain
  %----------------------------------------------------
  
  keep_going = 1;
  
  while keep_going
    
    % Redraw the image, BFC and XNot domain
    figure(1);
    imagesc(s); axis equal off; title(sprintf('%s : %d/%d',fname,ic,nims));
    hold on;
    plot(xbfc,ybfc,'og','markersize',12,'linewidth',3);
    circle(x0,y0,a0,'g','linewidth',3);
    plot(xx,yy,'r','linewidth',3);
    plot(pnts(1,:),pnts(2,:),'or','markersize',12,'linewidth',3);
    hold off;
    
    set(gcf,'Name','Define the long axis of XNot');
    
    % Initialize point vectors
    lax = zeros(2,1);
    lay = zeros(2,1);
    
    for p = 1:2
      
      % Get (x,y) coordinate from user
      [lax(p),lay(p),btn] = ginput(1);
      
      hold on;
      plot(lax,lay,'ow','markersize',12,'linewidth',3);
      hold off;
      
    end
    
    % Calculate great circle track passing through these two points
    [la,lo] = xy2sph(lax,lay,x0,y0,a0);
    [lagc,logc] = track2('gc',la(1),lo(1),la(2),lo(2));
    [xgc,ygc] = sph2xy(lagc,logc,x0,y0,a0);
    
    % Draw the XNot axis
    hold on;
    line(xgc,ygc,'color','w','linewidth',3);
    hold off;
    
    % Ask user to confirm XNot axis
    ans = questdlg('Confirmation','Accept XNot axis?','Yes','No','Yes');
    
    switch ans
      case 'Yes'
        keep_going = 0;
      case 'No'
        keep_going = 1;
    end
    
  end

  %------------------------------------------------------------
  % XNot DOMAIN ANALYSIS IN SPHERICAL COORDINATES
  %------------------------------------------------------------
  res = xnotanal(xx,yy,x0,y0,a0,lax,lay);

  fprintf('------------------------------------------\n');
  fprintf('Filename : %s\n',fname);
  fprintf('Image    : %d\n', ic);
  fprintf('------------------------------------------\n');
  % fprintf('XNot:\n');
  % fprintf('  Projected area : %0.3f\n', res.xnot_area);
  % fprintf('  Surface area   : %0.3f\n', res.xnot_surfarea);
  % fprintf('Embryo:\n');
  % fprintf('  Projected area : %0.3f\n', res.emb_area);
  % fprintf('  Surface area   : %0.3f\n', res.emb_surfarea);
  % fprintf('------------------------------------------\n');
  % fprintf('Projected Area Ratio : %0.3f\n', res.xnot_area / res.emb_area);
  fprintf('Surface Area Ratio   : %0.3f\n', res.xnot_surfarea / res.emb_surfarea);
  fprintf('------------------------------------------\n');
  fprintf('XNot length          : %0.3f deg\n', res.x_len);
  % fprintf('XNot min width       : %0.3f deg\n', res.y_min);
  fprintf('XNot max width       : %0.3f deg\n', res.y_max);
  fprintf('XNot mean width      : %0.3f deg\n', res.y_mean);
  fprintf('XNot median width    : %0.3f deg\n', res.y_med);
  fprintf('------------------------------------------\n');
  fprintf('XNot len/med(width)  : %0.3f\n', res.x_y);
  fprintf('\n');
  
  % Append results to file
  fprintf(fd, '%d ', ic);
  fprintf(fd, '%0.3f ', res.xnot_surfarea / res.emb_surfarea);
  fprintf(fd, '%0.3f ', res.x_y);
  fprintf(fd, '%0.3f ', res.xnot_surfarea);
  fprintf(fd, '%0.3f ', res.emb_surfarea);
  fprintf(fd, '%0.3f ', res.xnot_area);
  fprintf(fd, '%0.3f ', res.emb_area);
  fprintf(fd, '%0.3f ', res.x_len);
  fprintf(fd, '%0.3f ', res.y_min);
  fprintf(fd, '%0.3f ', res.y_max);
  fprintf(fd, '%0.3f ', res.y_mean);
  fprintf(fd, '%0.3f ', res.y_med);
  fprintf(fd, '%0.3f %0.3f %0.3f ', x0, y0, a0);
  fprintf(fd, '%0.3f ', pnts(1,:)); % x points for XNot 
  fprintf(fd, '%0.3f ', pnts(2,:)); % y points for XNot
  fprintf(fd, '%0.3f ', lax); % XNot axis lats (A,B)
  fprintf(fd, '%0.3f ', lay); % XNot axis lons (A,B)
  fprintf(fd, '\n');

end % Image loop

% Close results file
fclose(fd);