%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FilmFlow.m: Planar film low test case
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global b;       % number of lattice dirs
global weight;  % lattice weights
global cx;      % stencil x
global cy;      % stencil y
global inv;     % inverse directions
global lx;      % domain size x
global ly;      % domain size y


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   7   3   6 
%    \  |  /   
%   4 - 1 - 2 
%    /  |  \    
%   8   5   9 
%
% D2Q9 lattice definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b = 9;
weight = [4/9, 1/9,1/9,1/9,1/9, 1/36,1/36,1/36,1/36]; 
cx  = [ 0, 1, 0, -1, 0, 1, -1, -1, 1]; 
cy  = [ 0, 0, 1, 0, -1, 1,  1, -1, -1];
%dir = [ 1, 2, 3, 4,  5, 6,  7,  8,  9];
inv = [ 1, 4, 5, 2,  3, 8,  9,  6,  7];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization starts here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lx = 10; ly = 8;
tau = 0.95;
omega = 1/tau;
fx = +0.0001;
fy = -0.000000;

lambda_even = -1/tau
lambda_odd = -1/tau;
Lambda = 3.0/16.0;
lambda_odd = ( 2.0*lambda_even + 4.0 )/( (4*Lambda - 1.0)*lambda_even - 2.0 )


% create array of PDFs
f = zeros(ly,lx,9);
for i=1:9
    f(:,:,i) = weight(i);
end
%f(10,10,i) = 1.5*f(10,10,i);
ux = 0.0*ones(ly,lx);
uy = 0.0*ones(ly,lx);
f = polynomialEq(ux,uy);

interface = zeros(ly,lx);
interface(ly,:) = 1.0;

% create array for the boundary (set to 1 for obstacle)
boundary = zeros(ly,lx);
boundary(1,:) = 1.0;
%boundary(ly,:) = 1.0;

% analytical profile for film flow
nu = (2*tau-1)/6
h = ly-2.0; % 2 buffer nodes
js = [1:1:h].-0.5 % "real world" coordinates
uIdeal = fx/nu*(h.-0.5*js).*js
%uIdeal = [0 uIdeal fx/nu*(h.-0.5*(h)).*(h) ]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:12500
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % evaluation and visualization
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if(mod(t-1,100)==0)
   [rho, ux, uy] = calcMacro(f);
    ux = ux + 0.5*fx;           % Correct momentum for evaluation
    uy = uy + 0.5*fy;
    rho(find(boundary)) = 1.0;
    rho(find(interface)) = 1.0;
    ux(find(boundary))  = 0.0;
    uy(find(boundary))  = 0.0;
    ux(find(interface))= 0.0;
    uy(find(interface))=0.0;
    
    figure(1);
    subplot(2,1,1);
    imagesc((rho-1)/3); %colorbar;
    set(gca,'YDir','normal');
    title(['pressure']);

    subplot(2,1,2);
    plot(2:ly, rho(2:ly,floor(lx/2)));
    %rho
    
    figure(2);
    subplot(1,2,1);

    %[0 uIdeal fx/nu*(h.-0.5*(h)).*(h) ]

    plot( [1.5 js.+1.5 ly-0.5], [0 uIdeal fx/nu*(h.-0.5*(h)).*(h) ], ...
        (1:ly-2).+1.0, ux(2:ly-1,floor(lx/2))', '*', ...
        (1:ly-2).+1.0, (ux(2:ly-1,floor(lx/2))'-uIdeal), '+');
    title(['velocity profile after ',num2str(t-1),' time steps']);
    legend('analytical profile', 'simulation', 'absolute error');
    
    % profile of absolute error
    error=(ux(2:ly-1,floor(lx/2))'-uIdeal)'
    %ux(2:ly-1,floor(lx/2))'
    
    subplot(1,2,2);
    quiver(ux,uy);
    
    ## figure(3);
    ## imagesc(mass);  colorbar;
    ## set(gca,'YDir','normal');
    ## title(['mass']);
    
    %drawnow;
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % lattice Boltzmann
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [rho, ux, uy] = calcMacro(f); % TODO: optimize, collision sweep does this anew but we need it here for interfaceBeforeStream
   [f, neq] = collideTRT(f, lambda_even, lambda_odd);
   %f = collideBGK(f, omega);
   %neq = f - project(f);
   f = externalForce(rho,f,fx,fy);
   f = boundariesBeforeStream(f, boundary);
   
   %neq = f - project(f);
   %neq1-neq2

   %f = interfaceBeforeStream(f, interface, ux, uy);
   %f = interfaceLiBeforeStream(f, neq, interface, 0.5, lambda_even, ones(ly,lx), ux, uy);
   f = pabBeforeStream(f, neq, interface, lambda_even, ones(ly,lx), ux, uy);
   
   fNew = stream(f);
   %gradMass = gradient(mass);
   f = fNew;

end
