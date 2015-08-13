%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% collideTRT.m
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

function [ fOut, NonEq ] = collideTRT( fIn, lambda_even, lambda_odd )
global inv b;
%COLLIDE Relaxation according to Ginzburgs TRT model
%   lambda_even controls viscosity
%   lambda_odd should be chosen as
%     -( 2.0*lambda_even + 4.0 )/( (4*Lambda - 1.0)*lambda_even - 2.0 )
%   where Lambda = 1/4 controls the magic
    NonEq = fIn - project(fIn);
    
%     for i=1:b
%         fEven(:,:,i) = 0.5*(fIn(:,:,i) + fIn(:,:,inv(i)));
%         fOdd(:,:,i) = 0.5*(fIn(:,:,i) - fIn(:,:,inv(i)));
%     end
%     
%     for i=1:b
%         fEqEven(:,:,i) = 0.5*(fEq(:,:,i) + fEq(:,:,inv(i)));
%         fEqOdd(:,:,i) = 0.5*(fEq(:,:,i) - fEq(:,:,inv(i)));
%     end

    for i=1:b
        NEven(:,:,i) = 0.5*(NonEq(:,:,i) + NonEq(:,:,inv(i)));
    end
    
    fOut = fIn + lambda_even*(NEven) + lambda_odd*(NonEq-NEven);
end

