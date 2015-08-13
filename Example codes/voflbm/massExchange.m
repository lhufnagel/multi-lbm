%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% massExchange.m: 
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
%function [ deltaM ] = massExchange( mass, fOld, fNew )
function [deltaM]   = massExchange(liquid,interface,fOld,fNew)
global lx ly cx cy b inv;
%MASSEXCHANGE Given the old (post-collision) pdf-field fOld and the new
% (streamed) pdf-field fNew
%   Compute the mass exchange from the streaming
% ATTENTION: boundary conditions must be set already!!
    deltaM = zeros(ly,lx);
    %fluid = mass > 0.0;
    fluid = liquid+interface > 0.0;

    for i=1:b
        %deltaM(:,:) = deltaM(:,:) + circshift(mass, [cy(i), cx(i)]).*fNew(:,:,i) -
	%mass(:,:).*fOld(:,:,inv(i));

        deltaM(:,:) = deltaM(:,:) + circshift(fluid, [cy(i),cx(i)]).*(fNew(:,:,i) - fOld(:,:,inv(i)));
    end

    deltaM = deltaM.*fluid;
    deltaM = deltaM.*(liquid < 1.0);

    %fNew(:,:,i) - fOld(:,:,inv(i));
    
end


