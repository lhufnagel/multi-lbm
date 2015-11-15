clear;
% D2Q9 Lattice-Boltzmann/Level-Set Code for Multiphase flow
% According to doi:10.1016/j.camwa.2009.02.005
%
% Authors:
% Daniel Zint, daniel.zint@fau.de,
% Lorenz Hufnagel, lorenz.hufnagel@fau.de
%
% Used Libraries:
% -Toolbox of Level Set Methods, Copyright © 2007 by Ian M. Mitchell.
% http://www.cs.ubc.ca/~mitchell/ToolboxLS
% -Polygon intersection function, by Douglas M. Schwarz, BSD
% http://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections

%% Physical setup
t_end = 5; % [s]
x_len = 1; % [m]
y_len = 1; % [m]
rho_phys(1) = 1.1; % [kg/m^3] % Mass-density of Fluid 1,
rho_phys(2) = 1; % [kg/m^3] e.g. ~1000 for Water, 1.2 for Air. In [LU] Cell-density is varying around 1!
visc(1) = 1/6;% [m^2/s] kinematic viscosity of fluid 1
visc(2) = 1/6;% [m^2/s] kinematic viscosity of fluid 2
sigma = 1e-6;  % [kg/s^2] surface tension between the two phases, e.g. ~ 76E-3 between water and air

use_periodic_x = 1; % Use periodic boundaries on left/right domain border?
% -> Otherwise No-Slip/Bounce-Back
use_periodic_y = 1; % Use periodic boundaries on top/bottom domain border?
% -> Otherwise No-Slip/Bounce-Back, optionally moving Top-boundary
lidVel = 2; % [m/s] x-Velocity of the top Lid

testcase = 'line'; %possible values: 'line', 'circle'
% !!NOTE!! Initial Inteface is defined below (Line 87), once Level-Set Grid was initialized

%% LBM Setup
lbm_g=LBM_Grid;
lbm_g.dx=min(x_len,y_len)/20; % [m]
lbm_g.dt= lbm_g.dx^2;% [s]  %Diffusive scaling..

lbm_g.lidVel = lbm_g.dt/lbm_g.dx * [lidVel; 0];
lbm_g.omega(1) = 1/(3*visc(1)*lbm_g.dt/lbm_g.dx^2 + 1/2);
lbm_g.omega(2) = 1/(3*visc(2)*lbm_g.dt/lbm_g.dx^2 + 1/2);

lbm_it = 10; % Number of LBM-iterations until level-set update
% This number should be chosen according to the following:
% (lbm_it * lbm_g.dt) * max(interface-velocity) < lbm_g.dx
% Such, that the interface is not advected further than one cell per one LBM-LSM alternation

if (max(lbm_g.omega) > 1.7)
  disp('WARNING: omega > 1.7; interface handling maybe instable! Change dx, dt or viscosity');
end

lbm_g.nx = round(x_len/lbm_g.dx) + 2; % -> Ghost layer
lbm_g.ny = round(y_len/lbm_g.dx) + 2; % -> Ghost layer

%init LBM data
lbm_g.cells = ones(lbm_g.nx,lbm_g.ny,9);
for i=1:9
  lbm_g.cells(:,:,i) = lbm_g.weights(i);
end

%% Set up Level-Set library
% !!! ATTENTION !!!
% If Matlab warns "Undefined function or variable 'processGrid'",
% You have to adapt the path in addPathToKernel.m according to your OS.
run('./addPathToKernel');

% Create the Level-Set-grid.
lsm_g.dim = 2;
lsm_g.dx = lbm_g.dx;
lsm_g.min = [0.5*lsm_g.dx; 0.5*lsm_g.dx];

if(use_periodic_x || use_periodic_y)
  lsm_g.max = [x_len - 0.5*lsm_g.dx; y_len - 0.5*lsm_g.dx];
  lsm_g.bdry = @addGhostPeriodic;
else
  lsm_g.max = [x_len - 0.5*lsm_g.dx; y_len - 0.5*lsm_g.dx];
  lsm_g.bdry = @addGhostExtrapolate;
end
lsm_g = processGrid(lsm_g);

%---------------------------------------------------------------------------
% Create initial Interface conditions
%   Note that in the periodic BC case, these initial conditions will not
%   be continuous across the boundary
%   In practice that little detail is ignored
%   data: Implicit surface function, discretely given at xs
switch(testcase)
  case 'circle'
    % Celltype 1: Outer Fluid
    % Celltype 2: Inner Fluid
    center = [ 0.5*x_len; .5*y_len];
    radius = 0.25*min(x_len,y_len);
    
    data = zeros(size(lsm_g.xs{1}));
    data = data + (lsm_g.xs{1} - center(1)).^2 + (lsm_g.xs{2} - center(2)).^2;
    data = sqrt(data) - radius;
  case 'line'
    % Horizontal seperation line
    % Celltype 1: Upper fluid
    % Celltype 2: Lower fluid

    if (use_periodic_y)
      disp('Warning, periodic y boundaries are not meaningful!')
    end

    if (lidVel == 0)
      disp('Warning, zero lid velocity!')
    end

    if (rho_phys(1) ~= rho_phys(2))
      disp('Warning, different densities. Analytic solution does not account for that!')
    end
    
    seperator_y = 0.75*y_len;
    data = lsm_g.xs{2} - seperator_y;
  otherwise
    error('Testcase does not exist: %s', testcase);
end
%---------------------------------------------------------------------------

% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = zeros(lbm_g.nx,lbm_g.ny,2);
schemeData.grid = lsm_g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy. (Level-Set)
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.

accuracy = 'medium';
switch(accuracy)
  case 'low'
    schemeData.derivFunc = @upwindFirstFirst;
    integratorFunc = @odeCFL1;
  case 'medium'
    schemeData.derivFunc = @upwindFirstENO2;
    integratorFunc = @odeCFL2;
  case 'high'
    schemeData.derivFunc = @upwindFirstENO3;
    integratorFunc = @odeCFL3;
  case 'veryHigh'
    schemeData.derivFunc = @upwindFirstWENO5;
    integratorFunc = @odeCFL3;
  otherwise
    error('Unknown accuracy level %s', accuracy);
end

% Initialize Level Set Display
f = figure(1);
hold off;

h = visualizeLevelSet(lsm_g, data, 'contour', 0);

hold on;
grid on;
axis(lsm_g.axis);
daspect([ 1 1 1 ]);

tNow = 0;

celltype_old = dataToCelltype(data);
celltype_old = [celltype_old(:,1), celltype_old, celltype_old(:,end)];
celltype_old = [celltype_old(end,:); celltype_old; celltype_old(1,:)];

err_vek = [];
mass_vek = [];
pressure_vek =[(rho_phys(2)-rho_phys(1))/3];

startTime = cputime;

%---------------------------------------------------------------------------

%% Main Time Loop
%[0, t_end) (subject to a little roundoff).
while(t_end - tNow > 100 * eps * t_end)
  
  %% Pass Data from Level Set to LBM
  
  % Contour line of the (data == zero) level set
  C=contourc(lsm_g.xs{2}(1,:),lsm_g.xs{1}(:,1),  data, [0 0]);
  C=C(:,2:end);

  % Get curvature
  [curvature, ~] = curvatureSecond(lsm_g, data);
  
  % get first order derivative (=> normals)
  deriv = zeros(lbm_g.nx-2,lbm_g.ny-2,2);
  [derivL,derivR] = upwindFirstENO3(lsm_g,data,1);
  deriv(:,:,1) = 0.5 * (derivL + derivR);
  [derivL,derivR] = upwindFirstENO3(lsm_g,data,2);
  deriv(:,:,2) = 0.5 * (derivL + derivR);
  
  % "Adding" Ghost layers here!
  % -> Easier access with LBM indices
  data = [data(:,1), data, data(:,end)];
  data = [data(end,:); data; data(1,:)];
  
  deriv = [deriv(:,1,:), deriv, deriv(:,end,:)];
  deriv = [deriv(end,:,:); deriv; deriv(1,:,:)];
  
  curvature = [curvature(:,1), curvature, curvature(:,end)];
  curvature = [curvature(end,:); curvature; curvature(1,:)];
  
  % LBM Refill-Algorithm for cells that changed their (fluid) type
  celltype = dataToCelltype(data);
  changed_type = celltype - celltype_old;
  [x, y] = find(changed_type);
  
  for i=1:length(x)
    if x(i) < 2 || x(i) > lbm_g.nx-1 || y(i) < 2 || y(i) > lbm_g.ny-1
      continue;
    end
    
    normal = [deriv(x(i),y(i),1);deriv(x(i),y(i),2)];
    normal = normal/norm(normal); % normalize n
    normal = -1*changed_type(x(i),y(i))*normal; % make it point inward
    
    % Find lbm-link with smallest angle to interface normal
    [~, i_c_imax] = max((normal(1).*lbm_g.c(1,:) + normal(2).*lbm_g.c(2,:))./lbm_g.c_len);
    c_imax = lbm_g.c(:,i_c_imax);
    
    if celltype(x(i),y(i)) ~= celltype(x(i) + c_imax(1),y(i) + c_imax(2))
      disp('Warning in refill Algorithm, fluid only single Grid cell thick');
      
      % What to do if both/no links are pointing into same celltype?!?!
      %return;
      %continue;
    end
    
    % Here, we found the link pointing inward, i.e. away from the interface
    % Now refill, i.e.
    % 1. Interpolate density and velocity in current cell from neighbours along this link direction
    
    %Avoid the presence of three Ghost layers, hence modulo-magic
    x1 = mod([x(i);y(i)] + c_imax - 2, [lbm_g.nx-2; lbm_g.ny-2]) + 2;
    x2 = mod([x(i);y(i)]+2*c_imax - 2, [lbm_g.nx-2; lbm_g.ny-2]) + 2;
    x3 = mod([x(i);y(i)]+3*c_imax - 2, [lbm_g.nx-2; lbm_g.ny-2]) + 2;
    
    % "ensure" that all three nodes, we interpolate from are REALLY interior
    % (Could fail for small bubbles, with radius < 3 grid-cells....)
    if (celltype(x2) ~= celltype(x1))
      x2 = x1;
    end
    if (celltype(x3) ~= celltype(x2))
      x3 = x2;
    end
    
    % rho and vel are up-to-date from previous LBM (before Level-Set)
    rho_interp = 3 * rho(x1(1),x1(2)) - ...
      3 * rho(x2(1), x2(2)) + rho(x3(1), x3(2));
    
    vel_interp = 3 * vel(x1(1), x1(2),:) - ...
      3 * vel(x2(1), x2(2),:) + vel(x3(1),x3(2),:);
    
    % 2. Calculate equilibrium with rho_interp and vel_interp
    equil=zeros(1,9);
    non_eq=zeros(1,9);
    
    for j=1:9
      cTimesU = lbm_g.c(1,j) * vel_interp(1) + lbm_g.c(2,j) * vel_interp(2);
      equil(j) =  lbm_g.weights(j) .* (rho_interp + ...
        3 .* cTimesU + ...
        9/2 .* (cTimesU).^2 - ...
        3/2 .* (vel_interp(1).^2 + vel_interp(2).^2));
    end
    
    % 3. Take non-equlibrium from nearest neightbour, i.e. the one Link c_imax is pointing to
    for j = 1:9
      cTimesU = lbm_g.c(1,j) * vel(x(i) + c_imax(1), y(i) + c_imax(2),1) +...
        lbm_g.c(2,j) * vel(x(i) + c_imax(1), y(i) + c_imax(2),2);
      
      non_eq(j) = lbm_g.cells(x(i) + c_imax(1), y(i) + c_imax(2), j) - ...
        lbm_g.weights(j) .* (rho(x(i) + c_imax(1), y(i) + c_imax(2)) + ...
        3 .* cTimesU + ...
        9/2 .* (cTimesU).^2 - ...
        3/2 .* (vel(x(i) + c_imax(1), y(i) + c_imax(2),1).^2 + vel(x(i) + c_imax(1), y(i) + c_imax(2),2).^2));
    end
    
    % 4. Reinitialize cell
    lbm_g.cells(x(i),y(i),:) = equil(:) + non_eq(:);
  end
  
  %% LBM loop
  for t=1:lbm_it
    % LBM collide
    % and calculation of f_neq = f - f_eq
    rho = zeros(lbm_g.nx,lbm_g.ny);
    vel = zeros(lbm_g.nx,lbm_g.ny,2);
    f_eq = zeros(lbm_g.nx,lbm_g.ny,9);
    f_neq = zeros(lbm_g.nx,lbm_g.ny,9);
    
    for i=1:9
      rho(:,:) = rho(:,:) + lbm_g.cells(:,:,i);
      vel(:,:,1) = vel(:,:,1) + lbm_g.c(1,i) * lbm_g.cells(:,:,i);
      vel(:,:,2) = vel(:,:,2) + lbm_g.c(2,i) * lbm_g.cells(:,:,i);
    end
    
    if (max(abs(vel(:))) > 0.1)
      disp('Warning: LBM may become instable, vel > 0.1 [LU]');
    end
    
    for i=1:9
      cTimesU = lbm_g.c(1,i) * vel(:,:,1) + lbm_g.c(2,i) * vel(:,:,2);
      f_eq(:,:,i) = lbm_g.weights(i) .* (rho(:,:) + ...
        3 .* cTimesU(:,:) + ...
        9/2 .* (cTimesU(:,:)).^2 - ...
        3/2 .* (vel(:,:,1).^2 + vel(:,:,2).^2));
      
      f_neq(:,:,i) = lbm_g.cells(:,:,i) - f_eq(:,:,i);
      
      lbm_g.cells(:,:,i) = lbm_g.cells(:,:,i) - lbm_g.omega(celltype(:,:)) .* (lbm_g.cells(:,:,i) - f_eq(:,:,i));
    end
    
    %Interface treatment PRE-stream
    
    %Copy post-collision cells to cells_new, because we use it temporarily in the following
    lbm_g.cells_new(:,:,:) = lbm_g.cells(:,:,:);
    
    for x=1:lbm_g.nx
     for y=1:lbm_g.ny
       for k = 2:9
          
          x_b = x+lbm_g.c(1,k); % Indices of "boundary" point
          y_b = y+lbm_g.c(2,k); % in opposite fluid
          
          if min(x_b,y_b) < 1 || x_b > lbm_g.nx || y_b > lbm_g.ny
            continue
          end
          
          if celltype(x,y) == celltype(x_b,y_b)
            continue
          end
          
          % Approximate q linearly: the position of interface along the LBM-Link
          q = data(x_b, y_b)/(data(x_b, y_b) - data(x,y));
          
          if (x > 2 && x < lbm_g.nx-1 && y > 2 && y < lbm_g.ny-1)
            x_2 = [lsm_g.xs{1}(x-1,y-1);lsm_g.xs{2}(x-1,y-1)];
            x_1 = [lsm_g.xs{1}(x_b-1,y_b-1); lsm_g.xs{2}(x_b-1,y_b-1)];
            
            %Build more accurate intersection of link and contour-line. Mitchel-Library does not provide that..
            [x_int,y_int]=intersections([x_1(1),x_2(1)],[x_1(2),x_2(2)],C(2,:),C(1,:),false);
            if (~isempty(x_int))
              q = norm(x_1-[x_int(1);y_int(1)])/norm(x_1-x_2);
            end
          end
          
          % add_term1, First order interpolated velocity between phases to ensure continuity
          vel_int = q*[vel(x,y,1) ; vel(x,y,2)] + (1-q)*[vel(x_b, y_b, 1); vel(x_b,y_b,2)];
          add_term1 = 6*lbm_g.weights(k) * lbm_g.c(:,lbm_g.invDir(k))' * vel_int; % 6 f^*_i c_i u
          
          % S^(k)
          % Subscript according to Thoemmes-Paper!
          % 2: The point/cell under consideration
          % 1: The point, the current intersecting link points to
          
          S_2 = zeros(2);  % Shear rate tensors
          S_1 = zeros(2);
          for l = 1:9
            f_neq_2 = f_neq(x,y,l);
            f_neq_1 = f_neq(x_b,y_b,l);
            S_2 = S_2 + (lbm_g.c(:,l) * lbm_g.c(:,l)') * f_neq_2;
            S_1 = S_1 + (lbm_g.c(:,l) * lbm_g.c(:,l)') * f_neq_1;
          end
          S_2 = -1.5 * lbm_g.omega(celltype(x,y)) * S_2;
          S_1 = -1.5 * lbm_g.omega(celltype(x_b,y_b)) * S_1;
          
          
          % Lambda_i, Orthogonal projector
          Lambda_i = lbm_g.c(:,k)*lbm_g.c(:,k)' - lbm_g.c(:,k)'*lbm_g.c(:,k)/2*eye(2);
          
          % Lambda_i : [S]
          S_average = (S_2+S_1)*0.5;
          % Normal, tangent and curvature from LS-Toolbox

          %normal = [deriv(x,y,1);deriv(x,y,2)];
          normal = q*[deriv(x,y,1);deriv(x,y,2)]+(1-q)*[deriv(x_b,y_b,1);deriv(x_b,y_b,2)];
          normal = normal/norm(normal);        % normal n
          tangent = [-normal(2);normal(1)];    % tangent t
          %kappa = curvature(x,y); % curvature
          kappa = q*curvature(x,y)+(1-q)*curvature(x_b,y_b); % curvature
          
          % dynamic viscosity -> mu = mass_dens * nu
          mu_2 = lbm_g.dt/lbm_g.dx^2*visc(celltype(x,y)) * rho_phys(celltype(x,y));
          mu_1 = lbm_g.dt/lbm_g.dx^2*visc(celltype(x_b,y_b)) * rho_phys(celltype(x_b,y_b));
          mu_average = (mu_2 + mu_1)*0.5;
          
          % ! Formula given with opposite sign in doi:10.1016/j.jcp.2008.10.032!
          p_jump = 1/3 *(rho(x,y)*rho_phys(celltype(x,y)) - rho(x_b,y_b) * rho_phys(celltype(x_b,y_b)));
          mu_jump = (mu_2 - mu_1);
          
          % Intermedia results
          S_jump_n_n = 1/(2*mu_average) * (p_jump - 2*sigma*kappa) - mu_jump/mu_average * trace(S_average * (normal*normal'));
          S_jump_n_t = -mu_jump/mu_average * trace(S_average * (normal*tangent')');
          
          Lambda_times_S_jump = S_jump_n_n * ((normal'*lbm_g.c(:,k))^2 - (lbm_g.c(:,k)'*lbm_g.c(:,k))/2) + ...
            2*S_jump_n_t*(normal'*lbm_g.c(:,k))*(tangent'*lbm_g.c(:,k));
          
          % Lambda_i : S^(2)
          Lambda_times_S_2 = trace(Lambda_i*S_2');
          
          % add_term2 = R_i
          Lambda_times_A = (-q)*(1-q)*Lambda_times_S_jump - (q-0.5)*Lambda_times_S_2;
          add_term2 = 6*lbm_g.weights(k)*Lambda_times_A;   % R_i
          
          lbm_g.cells_new(x_b,y_b, lbm_g.invDir(k)) = lbm_g.cells(x,y,k) + add_term1 + add_term2;
        end
      end
    end
    
    %Copy back the (PRE-Stream) interface-treated links from cells_new to cells
    lbm_g.cells(:,:,:) = lbm_g.cells_new(:,:,:);
    
    % LBM stream
    for i=1:9
      lbm_g.cells_new(:,:,i) = circshift(lbm_g.cells(:,:,i), [lbm_g.c(1,i),lbm_g.c(2,i),0]);
    end
    
    % Domain boundary handling (This is POST-stream)
    % North & South
    if (use_periodic_y)
      lbm_g.cells_new(1:lbm_g.nx, 2, 2) = lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny, 2);
      lbm_g.cells_new(1:lbm_g.nx, 2, 3) = lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny, 3);
      lbm_g.cells_new(1:lbm_g.nx, 2, 9) = lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny, 9);
      
      lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny-1,6) = lbm_g.cells_new(1:lbm_g.nx,1,6);
      lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny-1,5) = lbm_g.cells_new(1:lbm_g.nx,1,5);
      lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny-1,7) = lbm_g.cells_new(1:lbm_g.nx,1,7);
      
      %Guarantee periodicity
      lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny, :)     = lbm_g.cells_new(1:lbm_g.nx, 2, :);
      lbm_g.cells_new(1:lbm_g.nx, 1, :)     = lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny-1, :);
    else
      %Solid boundaries
      lbm_g.cells_new(1:lbm_g.nx, 2, 2)    = lbm_g.cells_new(1:lbm_g.nx,1,6);
      lbm_g.cells_new(2:lbm_g.nx, 2, 3)    = lbm_g.cells_new(1:(lbm_g.nx-1),1,7);
      lbm_g.cells_new(1:(lbm_g.nx-1), 2, 9)= lbm_g.cells_new(2:lbm_g.nx,1,5);
      
      %Moving no-slip boundaries
      lbm_g.cells_new(1:lbm_g.nx, lbm_g.ny-1, 6)    = lbm_g.cells_new(1:lbm_g.nx,lbm_g.ny,2) + 6 * lbm_g.weights(2) * lbm_g.c(:,6)'*lbm_g.lidVel;
      lbm_g.cells_new(2:lbm_g.nx,  lbm_g.ny-1, 5)   = lbm_g.cells_new(1:(lbm_g.nx-1),lbm_g.ny,9) + 6 * lbm_g.weights(9) * lbm_g.c(:,5)'*lbm_g.lidVel;
      lbm_g.cells_new(1:(lbm_g.nx-1), lbm_g.ny-1, 7)= lbm_g.cells_new(2:lbm_g.nx,lbm_g.ny,3) + 6 * lbm_g.weights(3) * lbm_g.c(:,7)'*lbm_g.lidVel;
    end
    
    % East & west
    if (use_periodic_x)
      %periodic boundaries
      lbm_g.cells_new(2, 2:(lbm_g.ny-1), 4)  = lbm_g.cells_new(lbm_g.nx, 2:(lbm_g.ny-1), 4);
      lbm_g.cells_new(2, 3:(lbm_g.ny-1), 3)  = lbm_g.cells_new(lbm_g.nx, 3:(lbm_g.ny-1), 3);
      lbm_g.cells_new(2, 2:(lbm_g.ny-2), 5)  = lbm_g.cells_new(lbm_g.nx, 2:(lbm_g.ny-2), 5);

      lbm_g.cells_new(lbm_g.nx-1, 2:(lbm_g.ny-1), 8) = lbm_g.cells_new(1, 2:(lbm_g.ny-1), 8);
      lbm_g.cells_new(lbm_g.nx-1, 3:(lbm_g.ny-1), 9) = lbm_g.cells_new(1, 3:(lbm_g.ny-1), 9);
      lbm_g.cells_new(lbm_g.nx-1, 2:(lbm_g.ny-2), 7) = lbm_g.cells_new(1, 2:(lbm_g.ny-2), 7);
      
      %Guarantee periodicity
      lbm_g.cells_new(lbm_g.nx, 2:(lbm_g.ny-1), :)     = lbm_g.cells_new(2,2:(lbm_g.ny-1), :);
      lbm_g.cells_new(1, 2:(lbm_g.ny-1), :)  = lbm_g.cells_new(lbm_g.nx-1,2:(lbm_g.ny-1), :);
    else
      %Solid boundaries
      lbm_g.cells_new(lbm_g.nx-1, 2:(lbm_g.ny-1) ,8) = lbm_g.cells_new(lbm_g.nx,2:(lbm_g.ny-1), 4);
      lbm_g.cells_new(lbm_g.nx-1, 1:(lbm_g.ny-2), 7) = lbm_g.cells_new(lbm_g.nx,2:(lbm_g.ny-1), 3);
      lbm_g.cells_new(lbm_g.nx-1, 3:lbm_g.ny, 9) = lbm_g.cells_new(lbm_g.nx,2:(lbm_g.ny-1), 5);
      lbm_g.cells_new(2, 2:(lbm_g.ny-1)  ,4) = lbm_g.cells_new(1,2:(lbm_g.ny-1), 8);
      lbm_g.cells_new(2, 1:(lbm_g.ny-2),5) = lbm_g.cells_new(1,2:(lbm_g.ny-1), 9);
      lbm_g.cells_new(2, 3:lbm_g.ny,3)     = lbm_g.cells_new(1,2:(lbm_g.ny-1), 7);
    end
    
    %copy cells_new to cells, -> swap old and new buffer
    lbm_g.cells(:,:,:) = lbm_g.cells_new(:,:,:);
  end
  
  %% LBM visualisation
  
  % Collect data for visualisation
  rho = zeros(lbm_g.nx,lbm_g.ny);
  vel = zeros(lbm_g.nx,lbm_g.ny,2);
  for i=1:9
    rho(:,:) = rho(:,:) + lbm_g.cells(:,:,i);
    vel(:,:,1) = vel(:,:,1) + lbm_g.c(1,i) * lbm_g.cells(:,:,i);
    vel(:,:,2) = vel(:,:,2) + lbm_g.c(2,i) * lbm_g.cells(:,:,i);
  end
  p1 = 1/3*rho_phys(1)*mean(rho(celltype==1));
  p2 = 1/3*rho_phys(2)*mean(rho(celltype==2));
  pressure_vek = [pressure_vek, p2-p1];
  rho_sum = sum(sum(rho(2:lbm_g.nx-1,2:lbm_g.ny-1)));
  mass_vek = [mass_vek, rho_sum];
  
  %Velocity
  figure(2);
  subplot(1,2,1);
  %Scale vectors, such that, 0.05 velocity in [LU] := 1 vector length
  quiver(1/0.05*vel( 2:lbm_g.nx-1 , 2:lbm_g.ny-1 ,1)',...
    1/0.05*vel(2:lbm_g.nx-1 ,2:lbm_g.ny-1,2)','AutoScale','off');
  
  title(['Velocity, Range: [0, ' num2str(max(abs(vel(:)))) ']']);
  subplot(1,2,2);
  %Density
  contourf(rho(2:lbm_g.nx-1,2:lbm_g.ny-1)')
  title('Density');
  colorbar;
  
  % There are small mass losses. 
  % As mentioned in the literature, these stem from the Level-Set Method

  % figure(3);
  % plot(mass_vek/((lbm_g.nx-2)*(lbm_g.ny-2))-ones(size(mass_vek)));
  % title(['\Delta Mass (relative) over iterations']);
  
  if strcmp(testcase,'line')
    % Couette-Error
    
    a2 = lbm_g.lidVel(1)/(seperator_y + (visc(2)*rho_phys(2)/(visc(1)*rho_phys(1)))*(y_len - seperator_y));
    a1 = rho_phys(2)*visc(2)/(visc(1)*rho_phys(1)) * a2;
    offset = lbm_g.lidVel(1) - a1*y_len;
    
    analytic_profile = [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*ceil(seperator_y/lbm_g.dx)] -.5*lbm_g.dx),...
      (a1 * ([lbm_g.dx*ceil(1+seperator_y/lbm_g.dx+eps) : lbm_g.dx : y_len]-.5*lbm_g.dx)) + offset];
    
    lbm_profile = sqrt(vel(round(lbm_g.nx/2),2:lbm_g.ny-1,1).^2+vel(round(lbm_g.nx/2),2:lbm_g.ny-1,2).^2);
    
    figure(4)
    
    plot(lbm_profile,1:lbm_g.ny-2,'ro-');
    hold on
    
    plot(analytic_profile,1:lbm_g.ny-2,'b-');
    
    err_rel = norm((lbm_profile - analytic_profile)/analytic_profile)
    err_vek = [err_vek, err_rel];
    
    xlabel('Velocity [LU]');
    ylabel('y, Cell index');
    title('abs(velocity) at x-center column');
    legend({'LBM','analytic'},'Location','SouthEast')
    hold off
    
    figure(5)
    plot(err_vek);
    title(['(Relative) Error in L_2-Norm over \Delta T (:=' num2str(lbm_it) ' LBM-Iterations each)']);
  end
  
 if (strcmp(testcase,'circle'))
   % Pressure Drop
   figure(4)
   plot(pressure_vek+(rho_phys(1)-rho_phys(2))/3); %Offset by p_0
   hold on
   pressure_jump_analytic = 2*sigma/radius+(rho_phys(1)-rho_phys(2))/3;
   plot(pressure_jump_analytic*ones(size(pressure_vek)),'k--');
   title(['\Delta P over \Delta T (:=' num2str(lbm_it) ' LBM-Iterations each)']);
   %axis([0, length(pressure_vek), 0,2*pressure_jump_analytic]);
   legend({'LBM','analytic'},'Location','SouthEast')
   xlabel('10 LBM-Iterations');
   ylabel('\Delta P');
   hold off
   
 end
  
  celltype_old = celltype;
  %---------------------------------------------------------------------------
  
  %% Pass LBM data to Levelset
  % "remove" ghost layers
  schemeData.velocity = { lbm_g.dx/lbm_g.dt*vel(2:lbm_g.nx-1, 2:lbm_g.ny-1, 1); ...
    lbm_g.dx/lbm_g.dt*vel(2:lbm_g.nx-1, 2:lbm_g.ny-1, 2)};
  data = data(2:lbm_g.nx-1, 2:lbm_g.ny-1);
  
  %% Level Set code
  
  % Reshape data array into column vector for ode solver call.
  y0 = data(:);
  
  % How far to step?
  tSpan = [ tNow, min(t_end, tNow + lbm_g.dt * lbm_it) ];
  
  % Take a timestep.
  [ t, y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
    integratorOptions, schemeData);
  tNow = t(end);
  
  % Get back the correctly shaped data array
  data = reshape(y, lsm_g.shape);
  
  figure(f);
  [ figure_az, figure_el ] = view;
  
  delete(h);
  % Interface visualization.
  h = visualizeLevelSet(lsm_g, data, 'contour', 0, [ 't = ' num2str(tNow) ]);

  view(figure_az, figure_el);
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);

