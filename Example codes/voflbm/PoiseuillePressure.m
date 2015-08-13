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
lx = 8; ly = 10;
tau = 2.0;
omega = 1/tau;
fx = 0.0;

lambda_even = -1/tau
lambda_odd = -1/tau;
Lambda = 3.0/16.0;
lambda_odd = ( 2.0*lambda_even + 4.0 )/( (4*Lambda - 1.0)*lambda_even - 2.0 )


% create array of PDFs
f = zeros(ly,lx,9);
for i=1:9
    f(:,:,i) = weight(i);
end

rhoBig=1.1;
% create array for the boundary
boundary = zeros(ly,lx);
boundary(1,1:lx) = 1.0;
boundary(ly,1:lx) = 1.0;
pab = zeros(ly,lx);
pab(2:ly-1,1) = 1.0;
pab(2:ly-1,lx) = 1.0;
rhoPab = zeros(ly,lx);
rhoPab(2:ly-1,1) = rhoBig;
rhoPab(:,lx) = 1.0;

## boundary = zeros(ly,lx);
## boundary(1,:) = 1.0;
## boundary(ly,:) = 1.0;
## pab = zeros(ly,lx);
## pab(2:ly-1,1) = 1.0;
## pab(2:ly-1,lx) = 1.0;
## rhoPab = zeros(ly,lx);
## rhoPab(1:ly,1) = rhoBig;
## rhoPab(1:ly,lx) = 1.0;

% analytical profile
Force=fx+(rhoPab(2,1)-rhoPab(2,lx))/(3*(lx-2));
nu = (2*tau-1)/6;
n = ly-1;
L = ly-1; % subtract top and bottom wall
Uc = L*L*Force/(8*nu);
js = [1:1:n-1];
js-floor(ly/2);
uIdeal = Uc/((n)^2)*(2*js-1).*(2*n - 2*js - 1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:801
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % evaluation and visualization
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if(mod(t,100)==0)
    [rho, ux, uy] = calcMacro(f);
    ux = ux + 0.5*fx;           % Correct momentum for evaluation
    rho(find(boundary)) = 1.0;
    rho(find(pab)) = rhoPab(find(pab));
    ux(find(boundary))  = 0.0;
    uy(find(boundary))  = 0.0;
    ux(find(pab))  = 0.0;
    uy(find(pab))  = 0.0;
    
    figure(1);
    subplot(1,2,1);
    imagesc(rho); colorbar;

    subplot(1,2,2);
    plot([1.5, 2:1:lx-1 lx-0.5], [rhoBig rhoBig-3*Force*(0.5:lx-2)  1.0], '+-',
	[1:lx], rho(floor(ly/2),:), '*');
    title('pressure profile along center line');
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
    
    %uy(2,:)
    rho(floor(ly/2),:)
    
    subplot(1,2,2);
    quiver(ux,uy);
    drawnow;
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % lattice Boltzmann
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %f = collideBGK(f, omega);
   [rho, ux, uy] = calcMacro(f);
   [f, fNeq] = collideTRT(f, lambda_even, lambda_odd);
   f = externalForce(rho,f,fx,0);
   f = boundariesBeforeStream(f, boundary);
   f = pabBeforeStream(f, fNeq, pab, lambda_even, rhoPab, ux, uy);

   f = stream(f);
end
