%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% boundariesBeforeStream.m
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


function [ f ] = periodicBoundaries( f )
% periodic boundaries, BEFORE the streaming
 
   global b;    % number of lattice dirs
   global cx;   % stencil x
   global cy;   % stencil y
   global inv;  % inverse directions
   global lx;   % domain size x
   global ly;   % domain size y
   

   for k = 2:(ly-1)
       for dir = 1:b
           f(k,1,dir) = f(k,lx,dir); 
       end
   end
 
   for j = 1:lx
      for dir = 1:b
          f(1,j,dir)=f(ly,j,dir);
      end
   end

end
  
