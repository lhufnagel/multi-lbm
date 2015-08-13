%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% movingboundaries.m
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

function [ f ] = movingboundaries( f, boundary, rho, ub )
%BOUNDARIES does the half-way bounce back in cells
%   having a non-zero entry in the boundary field
%   rho is the density field
%   ub is a matrix prescribing the boundary velocity
   global b;        % number of lattice dirs
   global cx;       % stencil x
   global cy;       % stencil y
   global inv;      % inverse directions
   global weight;   % the lattice weights
   global lx;       % domain size x
   global ly;       % domain size y
   
   
   [j,i] = find(boundary);
   for k=1:length(i)
       for dir=1:b
           % let (ni,nj) be the boundary node with an off-domain node
           % in direction dir
           nj = mod(j(k)-cy(dir)-1,ly) + 1;
           ni = mod(i(k)-cx(dir)-1,lx) + 1;
           f(nj, ni, inv(dir)) = f(j(k), i(k), dir) ...
               + 6*weight(inv(dir))*rho(nj, ni)*(cx(inv(dir))*ub(j(k),i(k),1)+cy(inv(dir))*ub(j(k),i(k),2) );
           
       end
   end
end
