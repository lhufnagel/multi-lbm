%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pabBeforeStream.m
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

function [ f ] = interfaceLiBeforeStream( f, fNeq, boundary, delta, lambdaE, rho, ux, uy )
%Pressure-BOUNDARIES does the half-way pressure ANTI-bounce back according to Ginzburg in cells, BEFORE the streaming
% i.e. fetches the post-col pdfs from boundary nodes for the bc

   global b;    % number of lattice dirs
   global cx;   % stencil x
   global cy;   % stencil y
   global inv;  % inverse directions
   global lx;   % domain size x
   global ly;   % domain size y
   global weight; % lattice weights

   for i=1:b
        NEven(:,:,i) = 0.5*(fNeq(:,:,i) + fNeq(:,:,inv(i)));
	NOdd(:,:,i) = fNeq(:,:,i) - NEven(:,:,i);
   end

   [j,i] = find(boundary);
   for k=1:length(i)
       for dir=1:b
           % let (ni,nj) be the boundary node neighbor with the off-domain node
           % in direction dir
           nj = mod(j(k)-cy(dir)-1,ly) + 1;
           ni = mod(i(k)-cx(dir)-1,lx) + 1;
	   
	   nnj = mod(j(k)-2*cy(dir)-1,ly) + 1;
	   nni = mod(i(k)-2*cx(dir)-1,lx) + 1;
           
	   % interpolation coefficients
	   alpha = 1.0;
	   kappa1 = 1.0 - alpha*(0.5+delta);
	   kappam1 = 1.0 - 0.5*alpha;
	   kappa0 = alpha*delta -1.0;
	   C = alpha*(lambdaE*(0.5+delta)) - 2.0*lambdaE;

           % compute equilibrium+ term (only for non-linear equilibrium+
           %cu1 = 3*(cx(dir)*ux(nj,ni)+cy(dir)*uy(nj,ni));

	   % impose boundary condition for linear equilibrium
           f(j(k), i(k), inv(dir)) =  ...
	       alpha * weight(dir) * rho(j(k),i(k)) + ... % pressure condition
	       + C * NEven(nj,ni,dir) ... % shear stress condition (correction term)
	       + kappa1 * f(nj,ni,dir) ...
	       + kappam1 * f(nj,ni,inv(dir)) ...
	       + kappa0 * f(nnj,nni,dir);
       end
   end
end
