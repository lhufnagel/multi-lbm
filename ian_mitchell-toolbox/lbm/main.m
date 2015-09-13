function [ data, g, data0 ] = main(accuracy, displayType)
% convectionDemo: demonstrate a simple convective flow field.
%
%   [ data, g, data0 ] = convectionDemo(accuracy, displayType)
%  
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   displayType  String to specify how to display results.
%                  The specific string depends on the grid dimension;
%                  look at the helper visualizeLevelSet to see the options
%                  (optional, default depends on grid dimension).
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/9/04

t_end = 10; % [s]
x_len = 1; % [m]
y_len = 1; % [m]
lidVel = 2; % [m/s]
visc  = 1e-4;% [m^2/s]
lbm_it = 20; % No iterations until level-set update

lbm_g=Grid;
lbm_g.dx=0.02; % [m]
lbm_g.dt=0.001; % [s]
lbm_g.c_s=sqrt(1/3); % [m/s]
lbm_g.lidVel = [lidVel*lbm_g.dt/lbm_g.dx; 0];
lbm_g.omega = 1/(visc/(lbm_g.c_s^2*lbm_g.dt) + 1/2);
lbm_g.nx = x_len/lbm_g.dx + 2; %Ghost layer
lbm_g.ny = y_len/lbm_g.dx + 2; %Ghost layer

%init
lbm_g.cells = ones(lbm_g.nx,lbm_g.ny,9);
for i=1:9
  lbm_g.cells(:,:,i) = lbm_g.weights(i);
end

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('./addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = t_end;                  % End time.
plotSteps = 1/(lbm_g.dt*lbm_it);         % How many intermediate plots to produce?
t0 = 0;                      % Start time.
singleStep = 0;              % Plot at each timestep (overrides tPlot).

% Period at which intermediate plots should be produced.
tPlot = (tMax - t0) / (plotSteps - 1);

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
% What level set should we view?
level = 0;
% Pause after each plot?
pauseAfterPlot = 0;
% Delete previous plot before showing next?
deleteLastPlot = 0;
% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 1;

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
periodic = 1;

% Create the grid.
g.dim = 2;
g.min = 0;
g.dx = lbm_g.dx;
if(periodic)
  g.max = (1 - g.dx);
  g.bdry = @addGhostPeriodic;
else
  g.max = +1;
  g.bdry = @addGhostExtrapolate;
end
g = processGrid(g);

%---------------------------------------------------------------------------
% Choose the flow field.
v = @switchValue;

%---------------------------------------------------------------------------
% What kind of display?
if(nargin < 2)
  switch(g.dim)
   case 1
    displayType = 'plot';
   case 2
    displayType = 'contour';    
   case 3
    displayType = 'surface';
   otherwise
    error('Default display type undefined for dimension %d', g.dim);
  end
end

%---------------------------------------------------------------------------
% Create initial conditions (a circle/sphere).
%   Note that in the periodic BC case, these initial conditions will not
%   be continuous across the boundary unless the circle is perfectly centered.
%   In practice, we'll just ignore that little detail.
center = [ 0.45; .5; 0.0; 0.0 ];
radius = 0.35;
data = zeros(size(g.xs{1}));
for i = 1 : g.dim
  data = data + (g.xs{i} - center(i)).^2;
end
data = sqrt(data) - radius;
data0 = data;

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'low';
end

% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = v;
schemeData.grid = g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
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

if(singleStep)
  integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
end

%---------------------------------------------------------------------------
% Initialize Display
f = figure(1);

% Set up subplot parameters if necessary.
if(useSubplots)
  rows = ceil(sqrt(plotSteps));
  cols = ceil(plotSteps / rows);
  plotNum = 1;
  subplot(rows, cols, plotNum);
end

h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(t0) ]);

hold on;
if(g.dim > 1)
  axis(g.axis);
  daspect([ 1 1 1 ]);
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)
  %% Pass Data from level set to LBM
  celltype = dataToCelltype(data);
  %% LBM loop
  for t=1:lbm_it
    %lbm stream
    for i=1:9 
        lbm_g.cells_new(:,:,i) = circshift(lbm_g.cells(:,:,i), [lbm_g.c(1,i),lbm_g.c(2,i),0]); 
    end

    %lbm ugly boundary handling
    for x=1:lbm_g.nx
      %north & South
      lbm_g.cells_new(x  , 2,2) = lbm_g.cells_new(x,1,6);
      lbm_g.cells_new(x  , lbm_g.ny-1,6) = lbm_g.cells_new(x,lbm_g.ny,2) - 2/(lbm_g.c_s^2) * lbm_g.weights(2) * lbm_g.c(:,2)'*lbm_g.lidVel;
      if (x<lbm_g.nx)
        lbm_g.cells_new(x+1, 2,3) = lbm_g.cells_new(x,1,7);
        lbm_g.cells_new(x+1, lbm_g.ny-1,5) = lbm_g.cells_new(x,lbm_g.ny,9) - 2/(lbm_g.c_s^2) * lbm_g.weights(9) * lbm_g.c(:,9)'*lbm_g.lidVel;
      end
      if (x>1)
        lbm_g.cells_new(x-1, 2,9) = lbm_g.cells_new(x,1,5);
        lbm_g.cells_new(x-1, lbm_g.ny-1,7) = lbm_g.cells_new(x,lbm_g.ny,3) - 2/(lbm_g.c_s^2) * lbm_g.weights(3) * lbm_g.c(:,3)'*lbm_g.lidVel;
      end
    end

    for y=2:lbm_g.ny-1
      %east & west

      %solid boundaries
      %lbm_g.cells_new(lbm_g.nx-1, y ,8) = lbm_g.cells_new(lbm_g.nx,y, 4);
      %lbm_g.cells_new(lbm_g.nx-1, y-1, 7) = lbm_g.cells_new(lbm_g.nx,y, 3);
      %lbm_g.cells_new(lbm_g.nx-1, y+1, 9) = lbm_g.cells_new(lbm_g.nx,y, 5);
      %lbm_g.cells_new(2, y  ,4) = lbm_g.cells_new(1,y, 8);
      %lbm_g.cells_new(2, y-1,5) = lbm_g.cells_new(1,y, 9);
      %lbm_g.cells_new(2, y+1,3) = lbm_g.cells_new(1,y, 7);
      
      %periodic boundaries
      lbm_g.cells_new(2, y, 4)     = lbm_g.cells_new(lbm_g.nx,y, 4);
      lbm_g.cells_new(lbm_g.nx-1, y,8) = lbm_g.cells_new(1,y,8);
      if (y>2)
        lbm_g.cells_new(2, y, 3)     = lbm_g.cells_new(lbm_g.nx,y, 3);
        lbm_g.cells_new(lbm_g.nx-1, y,9) = lbm_g.cells_new(1,y,9);
      end
      if (y<lbm_g.ny-1)
        lbm_g.cells_new(2, y,5)     = lbm_g.cells_new(lbm_g.nx,y, 5);
        lbm_g.cells_new(lbm_g.nx-1, y,7) = lbm_g.cells_new(1,y,7);
      end
    end
    
    %boundary handling for the interface
    
    % calculation of diff_f = f - f_eq
    rho = zeros(lbm_g.nx,lbm_g.ny);
    vel = zeros(lbm_g.nx,lbm_g.ny,2);
    for i=1:9
      rho(:,:) = rho(:,:) + lbm_g.cells_new(:,:,i);
      vel(:,:,1) = vel(:,:,1) + lbm_g.c(1,i) * lbm_g.cells_new(:,:,i);
      vel(:,:,2) = vel(:,:,2) + lbm_g.c(2,i) * lbm_g.cells_new(:,:,i);
    end
    diff_f = zeros(lbm_g.nx,lbm_g.ny,9);
    for i = 1:9
        cTimesU = lbm_g.c(1,i) * vel(:,:,1) + lbm_g.c(2,i) * vel(:,:,2);
        diff_f(:,:,i) = lbm_g.cells_new(:,:,i) - ...
                    lbm_g.weights(i) .* (rho(:,:) + ...
                    1/(lbm_g.c_s^2) .* cTimesU(:,:) + ...
                    1/(2*lbm_g.c_s^4) .* (cTimesU(:,:)).^2 - ...
                    1/(2*lbm_g.c_s^2) .* (vel(:,:,1).^2 + vel(:,:,2).^2));
    end
    
    % get first order derivative
    deriv = zeros(lbm_g.nx-2,lbm_g.ny-2,2);
    [derivL,derivR] = upwindFirstENO3(g,data,1);
    deriv(:,:,1) = 0.5 * (derivL + derivR);
    [derivL,derivR] = upwindFirstENO3(g,data,2);
    deriv(:,:,2) = 0.5 * (derivL + derivR);
    
    [curvature, ~] = curvatureSecond(g, data);
    
    for x=1:lbm_g.nx-1-3
        for y = 1:lbm_g.ny-1-3
              for k = 2:9
                  if celltype(x+1,y+1) ~= celltype(x+1+lbm_g.c(1,k),y+1+lbm_g.c(2,k))
                    %% q
                    %q = 0.5;    % muss eigentlich aus level set berechnet werden
                    % pitfall: q und (1-q) müssen immer dem gleichen x zugeordnet werden
                    
                    q = data(x+1,y+1)/(data(x+1,y+1)-data(x+1+lbm_g.c(1,k),y+1+lbm_g.c(2,k))); 
                    % Die Interpolation müsste so stimmen, dadurch ist aber
                    % der add_term2 nicht mehr null, was das Programm zum
                    % Absturz bringt. Es ist zu vermuten, dass der Fehler
                    % in add_term2 steckt. Den muss ich eh nochmal
                    % aufräumen...
                    
                    
                    %% add_term1
                    vel_int = q*[vel(x,y,1) ; vel(x,y,2)] + (1-q)*[vel(x+lbm_g.c(1,k),y+lbm_g.c(2,k),1) ; vel(x+lbm_g.c(1,k),y+lbm_g.c(2,k),2)];
                    add_term1 = 6*lbm_g.dx*lbm_g.weights(k)*vel_int'*lbm_g.c(:,k); % 6 h f^*_i c_i
                    %% S^(k)
                    % Bei der Vergabe der Indizes habe ich mich an das Paper gehalten
                    % 2: Der Punkt den wir betrachten
                    % 1: Der Punkt im anderen Fluid
                    S_2 = zeros(2);   % S^(2)
                    S_1 = zeros(2);   % S^(1)
                    for l = 1:9
                        diff_f_2 = diff_f(x,y,l);
                        diff_f_1 = diff_f(x+lbm_g.c(1,k),y+lbm_g.c(2,k),l);
                        S_2 = S_2 + lbm_g.c(:,l) * lbm_g.c(:,l)' * diff_f_2;
                        S_1 = S_1 + lbm_g.c(:,l) * lbm_g.c(:,l)' * diff_f_1;
                    end
                    S_2 = -1.5 * lbm_g.omega * (1/lbm_g.dx^2) * S_2;
                    S_1 = -1.5 * lbm_g.omega * (1/lbm_g.dx^2) * S_1;
                    %% Lambda_i
                    Lambda_i = lbm_g.c(:,k)*lbm_g.c(:,k)' - (1.0/3.0)*norm(lbm_g.c(:,k))^2*eye(2);  % siehe S. 1143 oben
                    %% Lambda_i : [S] 
                    S_average = (S_2+S_1)*0.5;
                    % Normale, Tangente und Krümmung werden in der Toolbox bestimmt
                    normal = [deriv(x+1,y+1,1);deriv(x+1,y+1,2)];
                    normal = normal/norm(normal);        % normal n
                    tangent = [-normal(2);normal(1)];    % tangent t
                    kappa = curvature(x+1,y+1);          % curvature
                    
                    % mu = mass_dens * v
                    mu_2 = norm([vel(x,y,1) ; vel(x,y,2)]);
                    mu_1 = norm([vel(x+lbm_g.c(1,k),y+lbm_g.c(1,k),1) ; vel(x+lbm_g.c(1,k),y+lbm_g.c(1,k),2)]);
                    mu_average = (mu_2 + mu_1)*0.5;
                    mu_jump = mu_1 - mu_2;        % <-- kommt von mir. Leider konnte ich es nicht auf Richtigkeit prüfen
                    
                    p_jump = 1/(3*lbm_g.dx^2) *(rho(x+lbm_g.c(1,k),y+lbm_g.c(1,k)) - rho(x,y));    % Auf S. 1147 beschrieben
                    % In der Formel steckt noch die Massendichte. Diese
                    % habe ich vorerst außer acht gelassen, in dem
                    % Glauben, dass diese ungefähr 1 ist.
                    
                    %sigma = 0;          % surface tension
                    sigma = 0.016e-17;          % surface tension
                    
                    
                    % Zwischenergebnisse
                    S_jump_n_n = 1/(2*mu_average) * (p_jump + 2*sigma*kappa) - mu_jump/mu_average * trace(S_average * (normal*normal')');
                    S_jump_n_t = -mu_jump/mu_average * trace(S_average * (normal*tangent')');
                    
                    Lambda_times_S_jump = S_jump_n_n * ((normal'*lbm_g.c(:,k))^2 - (norm(lbm_g.c(:,k))^2)/3) + ...
                        2*S_jump_n_t*(normal'*lbm_g.c(:,k))*(tangent'*lbm_g.c(:,k));
                    
                    %% Lambda_i : S^(2)
                    Lambda_times_S_2 = trace(Lambda_i*S_2');
                    
                    %% add_term2 = R_i
                    Lambda_times_A = -q*(1-q)*Lambda_times_S_jump - (q-0.5)*Lambda_times_S_2;   % Teilergebnis zur Berechnung von R_i
                    add_term2 = 6*lbm_g.dx^2*lbm_g.weights(k)*Lambda_times_A;   % R_i
                    
                    %% do it
                    if k <= 5
                        lbm_g.cells_new(x+lbm_g.c(1,k),y+lbm_g.c(2,k),k+4) = lbm_g.cells_new(x,y,k) + add_term1 + add_term2;  
                    else
                        lbm_g.cells_new(x+lbm_g.c(1,k),y+lbm_g.c(2,k),k-4) = lbm_g.cells_new(x,y,k) + add_term1 + add_term2;  
                    end
                    
                  end
              end
        end
    end

    %lbm collide
    rho = zeros(lbm_g.nx,lbm_g.ny);
    vel = zeros(lbm_g.nx,lbm_g.ny,2);

    for i=1:9
      rho(:,:) = rho(:,:) + lbm_g.cells_new(:,:,i);
      vel(:,:,1) = vel(:,:,1) + lbm_g.c(1,i) * lbm_g.cells_new(:,:,i);
      vel(:,:,2) = vel(:,:,2) + lbm_g.c(2,i) * lbm_g.cells_new(:,:,i);
    end

    %! Here cells_new is implicitly swapped with the old cells -> Stream-Collide!
    for i=1:9
      cTimesU = lbm_g.c(1,i) * vel(:,:,1) + lbm_g.c(2,i) * vel(:,:,2);
      lbm_g.cells(:,:,i) = lbm_g.cells_new(:,:,i) - lbm_g.omega .* (lbm_g.cells_new(:,:,i) - ...
          lbm_g.weights(i) .* (rho(:,:) + ...
            1/(lbm_g.c_s^2) .* cTimesU(:,:) + ...
            1/(2*lbm_g.c_s^4) .* (cTimesU(:,:)).^2 - ...
            1/(2*lbm_g.c_s^2) .* (vel(:,:,1).^2 + vel(:,:,2).^2)));
    end
  end

  % LBM plot
  figure(2);
  % subplot(1,2,1);
  quiver(vel([2:lbm_g.nx-1],[2:lbm_g.ny-1],1)',...
   vel([2:lbm_g.nx-1],[2:lbm_g.ny-1],2)');
  % subplot(1,2,2);
  % contourf(rho([2:lbm_g.nx-1],[2:lbm_g.ny-1])')
  % colorbar;

  %% level set code
  
  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  if(pauseAfterPlot)
    % Wait for last plot to be digested.
    pause;
  end

  % Get correct figure, and remember its current view.
  figure(f);
  [ figure_az, figure_el ] = view;

  % Delete last visualization if necessary.
  if(deleteLastPlot)
    delete(h);
  end

  % Move to next subplot if necessary.
  if(useSubplots)
    plotNum = plotNum + 1;
    subplot(rows, cols, plotNum);
  end

  % Create new visualization.
  h = visualizeLevelSet(g, data, displayType, level, [ 't = ' num2str(tNow) ]);

  % Restore view.
  view(figure_az, figure_el);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);



%---------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------------------------
function out = switchValue(t, data, schemeData) %#ok<INUSL>
% switchValue: switches between two values.
%
%  out = switchValue(t, data, schemeData)
%
% Returns a constant value:
%           one     for t <= tSwitch;
%           two     for t >  tSwitch.
%
% By setting one and two correctly, this function can implement
%   the velocityFunc prototype for termConvection;
%   the scalarGridFunc prototype for termNormal, termCurvature and others;
%   and possibly some other prototypes...
%
% Parameters:
%   t            Current time.
%   data         Level set function.
%   schemeData   Structure (see below).
%
%   out          Either schemeData.one or schemeData.two.
%
% schemeData is a structure containing data specific to this type of 
%   term approximation.  For this function it contains the field(s)
%
%   .one         The value to return for t <= tSwitch.
%
% schemeData may contain other fields.

  %checkStructureFields(schemeData, 'one');

  out = {vel([2:lbm_g.nx-1],[2:lbm_g.ny-1],1); vel([2:lbm_g.nx-1],[2:lbm_g.ny-1],2)};
  end
end
