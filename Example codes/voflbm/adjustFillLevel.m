%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% adjustFillLevel.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample, Octave/Matlab script
% Copyright (C) 2014 Regina Ammer
% Address:  Chair for System Simulation,
%           Cauerstraße 11, 91058 Erlangen, Germany
% E-mail: regina.ammer@fau.de
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


function [ fill_level ] = adjustFillLevel( fill_level,mass,rho)
%BOUNDARIES does the half-way bounce back in cells, BEFORE the streaming
% i.e. copies the post-col pdfs from boundary nodes into the inverted
% direction of the obstacle nodes

   global b;    % number of lattice dirs
   global cx;   % stencil x
   global cy;   % stencil y
   global inv;  % inverse directions
   global lx;   % domain size x
   global ly;   % domain size y
   
   %[j,i] = find(interface);
   %for k=1:length(i)
    %   fill_level(i,j) = mass(i,j)./ rho(i,j);
   %end
   fill_level = mass./rho;
end
