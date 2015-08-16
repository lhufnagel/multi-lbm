%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% project.m
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

function fout = project(f)
% calculates the equilibrium for the whole lattice

global b

[rho,ux,uy] = calcMacro(f);
eqnorm = normeq(ux,uy);
fout = eqnorm;

for i=1:b
     fout(:,:,i)=rho.*eqnorm(:,:,i);
end