%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcMacro.m
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

function [ rho, ux, uy ] = calcMacro( f )
%CALCMACRO Summary of this function goes here
%   Detailed explanation goes here

    global b;
    global cx;
    global cy;
    
    rho = sum(f,3);
    ux  = (sum(f(:,:,[2 6 9]),3)-sum(f(:,:,[4 7 8]),3))./rho;
    uy  = (sum(f(:,:,[3 6 7]),3)-sum(f(:,:,[5 8 9]),3))./rho;
end

