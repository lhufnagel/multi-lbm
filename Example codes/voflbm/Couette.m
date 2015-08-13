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
lx = 16; ly = 14;
tau = 2.0;
omega = 1/tau;
fx = 0.0;

% create array of PDFs
f = zeros(ly,lx,9);
for i=1:9
    f(:,:,i) = weight(i);
end

% create array for the boundary
boundary = zeros(ly,lx);
ub = zeros(ly, lx, 2);

for i=1:lx
  for j=1:ly
      if((j-0.5) > (1+8 +0.25*(i-0.5)))
	boundary(j,i) = 1.0;
	ub(j,i,2) = 0.0025;
	ub(j,i,1) = 0.01;
      end
      if((j-0.5) < (1 + 0.25*(i-0.5)))
	boundary(j,i) = 1.0;
	ub(j,i,2) = -0.0025;
	ub(j,i,1) = -0.01;
      end
  end
end
    




% analytical profile
nu = (2*tau-1)/6;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:1000
    % Calculate Macroscopic
   [rho, ux, uy] = calcMacro(f);
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % evaluation and visualization
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if(t>1)
   if(mod(t,1)==0)
   	
    rho(find(boundary)) = 1.0;
    figure(1);
    imagesc(flipud(rho)); colorbar;
    
    figure(2);
    subplot(2,1,1);
    plot( ux(:,floor(lx/2))', '*');
    title(['velocity profile after ',num2str(t-1),' time steps']);
    legend( 'simulation');
    
    subplot(2,1,2);
    ux(find(boundary)) = 0.0;
    uy(find(boundary)) = 0.0;
    norm = ux.^2 + uy.^2;
    norm = sqrt(norm);
    quiver(ux./norm,uy./norm);
    drawnow;
   end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % lattice Boltzmann
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   f = collideBGK(f, omega);
   %f = externalForce(rho,f,fx,0);
   f = stream(f);
   f = periodicShift(f, -4);
   f = movingboundariesShifted(f, boundary, rho, ub, -4);
   

end
