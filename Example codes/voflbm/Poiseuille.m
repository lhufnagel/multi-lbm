%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Poiseuille.m: Poiseuille test case
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
lx = 8; ly = 13;
tau = 2.0;
omega = 1/tau;
fx = 0.0001;
fy = -0.0000000;

lambda_even = -1/tau
lambda_odd = -1/tau;
Lambda = 3.0/16.0;
lambda_odd = ( 2.0*lambda_even + 4.0 )/( (4*Lambda - 1.0)*lambda_even - 2.0 )


% create array of PDFs
f = zeros(ly,lx,9);
for i=1:9
    f(:,:,i) = weight(i);
end

% create array for the boundary
boundary = zeros(ly,lx);
boundary(1,:) = 1.0;
boundary(ly,:) = 1.0;
%boundary(10,10) = 1.0;

% analytical profile
nu = (2*tau-1)/6
n = ly-1;
L = ly-1; % subtract top and bottom wall
Uc = L*L*fx/(8*nu)
js = [1:1:n-1];
js-floor(ly/2);
uIdeal = Uc/((n)^2)*(2*js-1).*(2*n - 2*js - 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:1001
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % evaluation and visualization
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if(mod(t-1,100)==0)
    [rho, ux, uy] = calcMacro(f);
    ux = ux + 0.5*fx;           % Correct momentum for evaluation
    uy = uy + 0.5*fy;
    rho(find(boundary)) = 1.0;
    ux(find(boundary))  = 0.0;
    uy(find(boundary))  = 0.0;
    
    figure(1);
    %subplot(3,1,1);
    imagesc(rho); colorbar;
    %rho
    
    figure(2);
    subplot(1,2,1);
    
    plot( [0.5 js n-0.5]+1, [0 uIdeal 0], ...
        (1:ly-2)+1, ux(2:ly-1,floor(lx/2))', '*', ...
        (1:ly-2)+1, (ux(2:ly-1,floor(lx/2))'-uIdeal), '+');
    title(['velocity profile after ',num2str(t-1),' time steps']);
    legend('analytical profile', 'simulation', 'absolute error');
    
    % profile of absolute error
    %(ux(2:ly-1,floor(lx/2))'-uIdeal)

    %uy(floor(ly/2), floor(lx/2))

    ux(2:ly-1)

    subplot(1,2,2);
    quiver(ux,uy);
    drawnow;
  end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % lattice Boltzmann
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %f = collideBGK(f, omega);
   f = collideTRT(f, lambda_even, lambda_odd);
   f = externalForce(rho,f,fx,fy);
   f= boundariesBeforeStream(f,boundary);
   f = stream(f);
%   f = boundariesAfterStream(f, boundary);
   

end
