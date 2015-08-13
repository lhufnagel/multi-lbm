%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% entropicEq.m
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

function [fEq] = entropicEq( ux, uy )
%ENTROPICEQ Summary of this function goes here
%   Detailed explanation goes here

    global b;           % the number of dirs
    global dir;         % the lattice directions
    global weight;     % the lattice weights
    global cx;
    global cy;


    for dir=1:b
        sqtrmX = sqrt(1 + 3*ux.*ux);
        sqtrmY = sqrt(1 + 3*uy.*uy);
        fEq(:,:,dir)  = weight(dir) ...
                   .* (2 - sqtrmX) .* ((2*ux + sqtrmX)./(1-ux)).^cx(dir)...
                   .* (2 - sqtrmY) .* ((2*uy + sqtrmY)./(1-uy)).^cy(dir);
    end
end

