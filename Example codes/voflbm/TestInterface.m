%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TestInterface.m: test free surface (Koerner et al., 2005)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Lattice Boltzmann sample, Octave/Matlab script
% Copyright (C) 2014 Regina Ammer
% Address:  Chair for System Simulation,
%           Cauerstraße 11, 91058 Erlangen, Germany
% E-mail: regina.ammer@fau.de
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
lx = 10; ly = 10;
tau = 0.85;
omega = 1/tau;
fx = +0.0000;
fy = +0.0000; % test without force!
%%%%%%%%%%%%%%%%%%%

Lambda = 3.0/16.0;

% create array of PDFs
f = zeros(ly,lx,9);
for i=1:9
    f(:,:,i) = weight(i);
end


ux = 0.0*ones(ly,lx);
uy = 0.01*ones(ly,lx);

f = polynomialEq(ux,uy);

interface       = zeros(ly,lx);
interface(5,:)  = 1.0;
interface(10,:) = 1.0;


fill_level                  = zeros(ly,lx);
%fill_level(find(interface)) = 0.5;
fill_level(5,:)            = 0.95;  % 0.95;
fill_level(10,:)           = 0.08;  % 0.3;
fill_level(1:4,:)          = 1.0;
fill_level

mass 		      = zeros(ly,lx);
%mass(find(interface)) = 0.8;
mass(5,:)      = 0.95;
mass(10,:)     = 0.08;
mass(1:4,:)    = 1.0;

gas 	       = zeros(ly,lx);
gas(6:9,1:lx)  = 1.0;
 
liquid        = zeros(ly,lx);
liquid(1:4,:) = 1.0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:2000
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % evaluation and visualization
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if(mod(t-1,1)==0)
   [rho, ux, uy] = calcMacro(f);


    ux = ux + 0.5*fx;           % Correct momentum for evaluation
    uy = uy + 0.5*fy;
    
    rho(find(interface)) = 1.0;
    rho(find(gas))       = 0.0;    
        
    figure(1);
    subplot(2,1,1);
    imagesc((rho-1)/3); %colorbar;
    set(gca,'YDir','normal');
    title(['pressure']);

    subplot(2,1,2);
    plot(2:ly, rho(2:ly,floor(lx/2)));
    title(['rho']);
      
    figure(2);
    quiver(ux,uy);
    
    figure(3);
    imagesc(mass);  colorbar;
    ## set(gca,'YDir','normal');
    title(['mass']);
    

    figure(4);
    imagesc(gas); colorbar;
    title(['gas']);
    %drawnow;
   end

   t
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % lattice Boltzmann
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   [rho, ux, uy] = calcMacro(f); % TODO: optimize, collision sweep does this anew but we need it here for interfaceBeforeStream

   f = collideBGK(f, omega);
   %f = externalForce(rho,f,fx,fy);
   f = interfaceBeforeStream( f, gas, ux, uy);
   %f = interfaceBeforeStream( f, gas, ones(ly,lx), -0.01*ones(ly,lx));

   fNew = stream(f);
   %fNew = boundariesAfterStream(fNew,boundary);
   
   conv_I_L = zeros(ly,lx);
   conv_I_G = zeros(ly,lx);
   
   % before=rho
   [mass,fill_level,liquid,interface,gas,fNew]= cellConversion4(ly,lx,mass,fill_level,liquid,interface,gas,conv_I_L,conv_I_G,fNew,b,weight,cx,cy);
   % after=rho

   me = massExchange(liquid, interface,f,fNew);
   mass = mass + me; 
   %%% here do cell conversions!
   f = fNew; 

   % adjust new fill_level
   fill_level = adjustFillLevel(fill_level,mass,rho); 
   
   
   %if(mod(t-1,10) == 1)
   	t
   	mass
   	fill_level
   	rho
   	liquid;
   	interface;
   	gas;
   %end
   
end
