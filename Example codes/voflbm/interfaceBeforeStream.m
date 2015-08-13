%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interfaceBeforeStream.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample, Octave/Matlab script
% Copyright (C) 2014 Simon Bogner
% Address:  Chair for System Simulation,
%           Cauerstraße 11, 91058 Erlangen, Germany
% E-mail: Simon.Bogner@fau.de
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public 
% License along with this program; if not, write to the Free 
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
% Boston, MA  02110-1301, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ f ] = interfaceBeforeStream( f, boundary, ux, uy )
%INTERFACE-BOUNDARIES does the half-way pressure bounce back according to Koerner in cells, BEFORE the streaming
% i.e. fetches the post-col pdfs from boundary nodes for the bc

   global b;    % number of lattice dirs
   global cx;   % stencil x
   global cy;   % stencil y
   global inv;  % inverse directions
   global lx;   % domain size x
   global ly;   % domain size y
   global weight; % lattice weights
   
   [j,i] = find(boundary);
   for k=1:length(i)
       for dir=1:b
           % let (ni,nj) be the boundary node neighbor with the off-domain node
           % in direction dir
           nj = mod(j(k)-cy(dir)-1,ly) + 1;
           ni = mod(i(k)-cx(dir)-1,lx) + 1;
           
           % compute equilibrium
           cu1 = 3*(cx(dir)*ux(nj,ni)+cy(dir)*uy(nj,ni));
           f(j(k), i(k), inv(dir)) =  2.0 * weight(dir) * ...
                       ( 1 + 1/2*(cu1*cu1)  - 3/2*(ux(nj,ni)^2+uy(nj,ni)^2) ) ...
                       -  f(nj, ni, dir);
       end
   end
end
