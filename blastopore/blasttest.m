function blasttest
% blasttest
%
% Compare projected and calculated surface areas using FrogSpawn functions.
%
% AUTHOR: Mike Tyszka, Ph.D.
% PLACE : Caltech BIC
% DATES : 03/15/2004 JMT From scratch
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

x0 = 0;
y0 = 0;

av = 10:10:100
hs_err = zeros(size(av));

for ac = 1:length(av)
  
  a = av(ac);
  
  [xm,ym] = meshgrid((-a:a) + x0,(-a:a) + y0);

  r = sqrt((xm-x0).^2 + (ym-y0).^2);

  [xx,yy] = circle(x0,y0,a);

  carea = polyarea(xx,yy);
  hsarea = spharea(xx,yy,x0,y0,a);
  act_carea = pi * a * a;
  act_hsarea = 2 * pi * a * a;

  fprintf('Circle area            : %0.1f\n', carea);
  fprintf('Hemisphere area        : %0.1f\n', hsarea);
  fprintf('Actual circle area     : %0.1f\n', act_carea);
  fprintf('Actual hemisphere area : %0.1f\n', act_hsarea);
  
  hs_err(ac) = (hsarea - act_hsarea) / act_hsarea * 100;

end

plot(av,hs_err);
ylabel('Error (%)');
set(gca,'YLim',[-0.1 0.1]);