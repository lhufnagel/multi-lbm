%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% periodicShift.m
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

function [ out ] = periodicShift( pdfs, number )
%STREAM Stream the distribution functions
    
    global b;
    global dir; % the lattice directions (stencil)
    global cx;
    global cy;
    global lx;
    global ly;

    out = pdfs;
    buffer1(:,:) = pdfs(:,1,:);
    buffer2(:,:) = pdfs(:,lx,:);
    ## for dir=1:b
    ##   out(:,:,dir) = circshift(pdfs(:,:,dir), [number,0]);
    ##   out(:,:,dir) = circshift(pdfs(:,:,dir), [-number,0]);
    ## end
    
    buffer1 = circshift(buffer1, [number,0]);
    buffer2 = circshift(buffer2, [-number,0]);
    
    out(:,1,[2 6 9]) = buffer1(:,[2 6 9]);
    out(:,lx,[4 7 8]) = buffer2(:,[4 7 8]);
end
