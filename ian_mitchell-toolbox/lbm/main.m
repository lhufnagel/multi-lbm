clear;
nargin = 0;
%
% Code based on convectionDemo

%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

t_end = 100; % [s]
x_len = 1; % [m] Achtung! Scheint hardgecodet im Level-Set zu sein
y_len = 1; % [m]
rho_phys(1) = 1; % [kg/m^3] % e.g. 1000 for Water, 1,2 for Air. In Lattice-Boltzmann-Units Cell-density is varying around 1. 
            % To obtain physical value we multiply by physical density, see e.g. formula for pressure jump, page 1147
rho_phys(2) = 1; % [kg/m^3] % 
lidVel = .1; % [m/s]
visc(1)  = 5e-1;% [m^2/s] % Oben im Line-Testcase, muss kleinere Viskosität haben als die untere Schicht, sonst sinnlos
visc(2)  = 3e-1;% [m^2/s]
sigma = 1.6e-4; %0.016; % [N/m] surface tension. Material parameter between different fluids, e.g. ~ 76*10^-3 between water and air 
lbm_it = 50; % Number of iterations until level-set update
% Diese Zahl sollte mMn folgendermaßen beschränkt sein:
% Delta_T = lbm_it * lbm_g.dt ist das Zeitinterval, in dem sich Level-Set und LBM abwechseln.
% Die maximale Geschwindigkeit in der Lid-Driven-Cavity ist lidVel (wenn man starke Krümmungseffekte und daraus folgende große Oberflächenspannungen an kleinen Blasen ignoriert).
% Da ich möchte, dass sich das Interface MAXIMAL eine Zelle weit bewegt pro Delta_T (CFL-Bedingung, außerdem geht sonst der Refill-Algorithmus kaputt),
% muss lidVel * Delta_T = lidVel * lbm_it * lbm_g.dt < lbm_g.dx sein!

err_vek = [];

lbm_g=Grid;
lbm_g.dx=0.05/x_len; % [m]

lbm_g.dt=0.01;% [s] 
lbm_g.lidVel = lbm_g.dt/lbm_g.dx * [lidVel; 0];
lbm_g.omega(1) = 1/(3*visc(1)*lbm_g.dt/lbm_g.dx^2 + 1/2); 
lbm_g.omega(2) = 1/(3*visc(2)*lbm_g.dt/lbm_g.dx^2 + 1/2);

if (max(lbm_g.omega) > 1.7 || min(lbm_g.omega) < 0.6)  %%TODO abklären? bound von unten: arXiv:0812.3242v2 letzter Absatz 1/1.8 = 0.555
  disp('WARNING: omega out of range; interface handling maybe instable! Change dx, dt or viscosity');
end

if (max(abs(lbm_g.lidVel)) > 0.1)
      disp('Warning: Top-Lid too fast, vel > 0.1 [LU]');
end

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
plotSteps = t_end/(lbm_g.dt*lbm_it) + 1;         % How many intermediate plots to produce?
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
deleteLastPlot = 1;
% Plot in separate subplots (set deleteLastPlot = 0 in this case)?
useSubplots = 0;

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

testcase = 'circle';
seperator_y = 0.42;

switch(testcase)
    case 'circle'
        center = [ 0.45; .5; 0.0; 0.0 ];
        radius = 0.35;
        data = zeros(size(g.xs{1}));
        for i = 1 : g.dim
          data = data + (g.xs{i} - center(i)).^2;
        end
        data = sqrt(data) - radius;
    case 'line'
        % Horizontale Trennlinie
        seperator_y = 0.42;
        data = g.xs{2} - seperator_y;
    otherwise
        error('Testcase does not exist: %s', testcase);
end
data0 = data;

%---------------------------------------------------------------------------
if(nargin < 1)
  accuracy = 'medium';
end

% Set up spatial approximation scheme.
schemeFunc = @termConvection;
schemeData.velocity = zeros(lbm_g.nx,lbm_g.ny,2);
schemeData.grid = g;

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy. (Level-Set)
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
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

h = visualizeLevelSet(g, data, 'contour', level, [ 't = ' num2str(t0) ]);

hold on;
if(g.dim > 1)
  axis(g.axis);
  daspect([ 1 1 1 ]);
end

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;

celltype_old = dataToCelltype(data);

startTime = cputime;
while(tMax - tNow > small * tMax)

  % Pass Data from level set to LBM

  % "Adding" Ghost layers here! 
  celltype = dataToCelltype(data);

  % Refill-Algorithm for cells that changed their type
  changed_type = celltype - celltype_old;
  [x, y] = find(changed_type);

  % get first order derivative (=> normals)
  deriv = zeros(lbm_g.nx-2,lbm_g.ny-2,2);
  [derivL,derivR] = upwindFirstENO3(g,data,1);
  deriv(:,:,1) = 0.5 * (derivL + derivR);
  [derivL,derivR] = upwindFirstENO3(g,data,2);
  deriv(:,:,2) = 0.5 * (derivL + derivR);
    
  for i=1:length(x)
    if x(i) < 3 || x(i) > g.N(1) || y(i) < 3 || y(i) > g.N(2)
      continue;
    end

    normal = [deriv(x(i)-1,y(i)-1,1);deriv(x(i)-1,y(i)-1,2)];
    normal = normal/norm(normal); % normalize n
    normal = -1*changed_type(x(i),y(i))*normal; % make it point inward

    % Find lbm-link with smallest angle to interface normal
    % -> Maximum of inner-product of normed vectors
    [~, i_c_imax] = max((normal(1).*lbm_g.c(1,:) + normal(2).*lbm_g.c(2,:))./lbm_g.c_len);
    c_imax = lbm_g.c(:,i_c_imax);

    if celltype(x(i),y(i)) ~= celltype(x(i) + c_imax(1),y(i) + c_imax(2))
      %TODO, what if both/no links are pointing into same celltype?!?!
      disp('bad error in refill Algorithm');
      return;
    end
    % We found the link pointing inward, i.e. away from the interface
    % Now refill, i.e. 
    % 1. interpolate density and velocity in current cell from neighbours along this link direction

    % TODO ensure that all three nodes, we interpolate from are REALLY interior 
    % (Could fail for small bubbles, with radius < 3 grid-cells....)

    rho_interp = 3 * rho(x(i) + c_imax(1), y(i) + c_imax(2)) - ...
       3 * rho(x(i) + 2*c_imax(1), y(i) + 2*c_imax(2)) + ...
       rho(x(i) + 3*c_imax(1), y(i) + 3*c_imax(2)); 

    % TODO Evtl aus nächstem Timestep interpolieren?!
    % UND: Ich verstehe nicht, wieso die im Paper den vektor abziehen... er zeigt doch nach innen und ich interpoliere von inneren Punkten..

    vel_interp = 3 * vel(x(i) + c_imax(1), y(i) + c_imax(2),:) - ...
       3 * vel(x(i) + 2*c_imax(1), y(i) + 2*c_imax(2),:) + ...
       vel(x(i) + 3*c_imax(1), y(i) + 3*c_imax(2),:); 

    % TODO hier wirds sehr unschön: In der Interpolationsformel im Thömmes-Paper wird für u_Gamma die Geschwindigkeit v aus dem Level-Set (?!?! Siehe auch Paper Kapitel 5.3) verwendet. 
    % Diese Geschwindigkeit stammt - so wie ich es verstehe - aus der Fast-Marching/Constant extension velocity  methode. 
    % Diese Methode implementiert die Bibliothek von Mitchel allerdings nicht. (Im Gegensatz zu evtl http://ktchu.serendipityresearch.org/software/lsmlib/ )
    % -> ?? 
    % ... -> ich verwende die Formel von rho auch für die Geschwindigkeit. 
    % Die Formel ist auch zu finden in Lallemand & Luo: Lattice Boltzmann method for moving boundaries. J Comput Phys 184(2): 406- 421
    % und sollte laut dem Paper ähnliche Ergebnisse produzieren (siehe Ende von Kapitel 3 in dem Paper)



    % 2. equilibrium berechnen mit rho_interp und vel_interp
    for j=1:9
      cTimesU = lbm_g.c(1,j) * vel_interp(1) + lbm_g.c(2,j) * vel_interp(2);
      equil(j) =  lbm_g.weights(j) .* (rho_interp + ...
            3 .* cTimesU + ...
            9/2 .* (cTimesU).^2 - ...
            3/2 .* (vel_interp(1).^2 + vel_interp(2).^2));
    end

    
    % 3. non_eq von nearest neightbour also zelle(x(i) + c_imax(1), y(i) + c_imax(2)) kopieren
    for j = 1:9
        cTimesU = lbm_g.c(1,j) * vel(x(i) + c_imax(1), y(i) + c_imax(2),1) + lbm_g.c(2,j) * vel(x(i) + c_imax(1), y(i) + c_imax(2),2);
        non_eq(j) = lbm_g.cells(x(i) + c_imax(1), y(i) + c_imax(2), j) - ...
                    lbm_g.weights(j) .* (rho(x(i) + c_imax(1), y(i) + c_imax(2)) + ...
                    3 .* cTimesU + ...
                    9/2 .* (cTimesU).^2 - ...
                    3/2 .* (vel(x(i) + c_imax(1), y(i) + c_imax(2),1).^2 + vel(x(i) + c_imax(1), y(i) + c_imax(2),2).^2));
    end

    % 4. Reinitialise
    lbm_g.cells(x(i),y(i),:) = equil(:) + non_eq(:);
  end
    

  rho_old = 0;

  %% LBM loop
  for t=1:lbm_it
    % lBM collide
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

    for i=1:9
      cTimesU = lbm_g.c(1,i) * vel(:,:,1) + lbm_g.c(2,i) * vel(:,:,2);
      f_eq(:,:,i) = lbm_g.weights(i) .* (rho(:,:) + ...
          3 .* cTimesU(:,:) + ...
          9/2 .* (cTimesU(:,:)).^2 - ...
          3/2 .* (vel(:,:,1).^2 + vel(:,:,2).^2));

      f_neq(:,:,i) = lbm_g.cells(:,:,i) - f_eq(:,:,i);
          
      lbm_g.cells(:,:,i) = lbm_g.cells(:,:,i) - lbm_g.omega(celltype(:,:)) .* (lbm_g.cells(:,:,i) - f_eq(:,:,i));
    end

    
    % boundary handling for the interface (PRE-Stream)
    
    % get first order derivative -> Surface normal
    deriv = zeros(lbm_g.nx-2,lbm_g.ny-2,2);
    [derivL,derivR] = upwindFirstENO3(g,data,1);
    deriv(:,:,1) = 0.5 * (derivL + derivR);
    [derivL,derivR] = upwindFirstENO3(g,data,2);
    deriv(:,:,2) = 0.5 * (derivL + derivR);
    
    % get curvature
    [curvature, ~] = curvatureSecond(g, data);
     
    lbm_g.cells_new(:,:,:) = lbm_g.cells(:,:,:);
    
    data = [data(:,1), data, data(:,end)];
    data = [data(end,:); data; data(1,:)];

    deriv = [deriv(:,1,:), deriv, deriv(:,end,:)];
    deriv = [deriv(end,:,:); deriv; deriv(1,:,:)];

    curvature = [curvature(:,1), curvature, curvature(:,end)];
    curvature = [curvature(end,:); curvature; curvature(1,:)];

    for x=1:lbm_g.nx
        for y = 1:lbm_g.ny
              for k = 2:9

                  if x+lbm_g.c(1,k) < 1 || x+lbm_g.c(1,k) > g.N(1)+2 || y+lbm_g.c(2,k) < 1 || y+lbm_g.c(2,k) > g.N(2)+2  %TODO evtl anpassen
                      continue
                  end

                  %Hier evtl die "isNearInterface(..)"-Funktion verwenden? (Doku S. 124)
                  if celltype(x,y) == celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))
                      continue
                  end

                  % q
                  q = data(x+lbm_g.c(1,k), y+lbm_g.c(2,k))/(data(x+lbm_g.c(1,k), y+lbm_g.c(2,k)) - data(x,y)); 
                  % Die Liniensegmente vom "exakten" Inteface kriegt man mit contourc(g.xs{1}(:,1), g.xs{2}(1,:), data, [0 0]);
                  % Rückgabe-Format Doku: http://de.mathworks.com/help/matlab/ref/contour-properties.html#prop_ContourMatrix
                  % Eventuell kann man damit q noch genauer berechnen (also exakt den Schnittpunkt mit dem Link c_i bilden)
                  % N.B.: Mit dem contourc-Aufruf arbeitet zumindest der visualizeLevelSet-Aus der Mitchel-Bibliothek intern

                  % Außerdem TODO : Vermutlich verliert data die  Signed-Distance-Eigenschaft. (Da wir nicht die Extension-Velocites generieren, i.e Fast-Marching-Methode nicht verwenden siehe Thömmes-Paper Kapitel 4.2.)
                  % Folglich brauchen wir einen Aufruf von signedDistanceIterative(..), siehe Kapitel 3.4.7 in ToolboxLS-Doku.
                  % Sonst wird hier q unsinnig berechnet(??)

                  %% add_term1
                  vel_int = q*[vel(x,y,1) ; vel(x,y,2)] + (1-q)*[vel(x + lbm_g.c(1,k), y+lbm_g.c(2,k), 1); vel(x  +lbm_g.c(1,k),y  +lbm_g.c(2,k),2)];
                  add_term1 = 6*lbm_g.weights(k) * lbm_g.c(:,k)' * vel_int; % 6 f^*_i c_i u

                  %% S^(k)
                  % Bei der Vergabe der Indizes habe ich mich an das Paper gehalten
                  % 2: Der Punkt den wir betrachten
                  % 1: Der Punkt im anderen Fluid
                  S_2 = zeros(2);   % S^(2)
                  S_1 = zeros(2);   % S^(1)
                  for l = 1:9
                      f_neq_2 = f_neq(x,y,l);
                      f_neq_1 = f_neq(x+lbm_g.c(1,k),y+lbm_g.c(2,k),l);
                     %S_2 = S_2 + (lbm_g.c(:,l) * lbm_g.c(:,l)' - lbm_g.c(:,l)'*lbm_g.c(:,l)/2*eye(2))* f_neq_2; %Enforce traceless-ness. Does not seem to have any impact
                     %S_1 = S_1 + (lbm_g.c(:,l) * lbm_g.c(:,l)' - lbm_g.c(:,l)'*lbm_g.c(:,l)/2*eye(2))* f_neq_1;
                      S_2 = S_2 + (lbm_g.c(:,l) * lbm_g.c(:,l)') * f_neq_2;
                      S_1 = S_1 + (lbm_g.c(:,l) * lbm_g.c(:,l)') * f_neq_1;
                  end
                  S_2 = -1.5 * lbm_g.omega(celltype(x,y)) * S_2;
                  S_1 = -1.5 * lbm_g.omega(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))) * S_1;

                  %% Lambda_i
                  Lambda_i = lbm_g.c(:,k)*lbm_g.c(:,k)' - lbm_g.c(:,k)'*lbm_g.c(:,k)/2*eye(2);  % siehe S. 1143 oben
                  %% Lambda_i : [S] 
                  S_average = (S_2+S_1)*0.5;
                  % Normale, Tangente und Krümmung werden in der Toolbox bestimmt
                  normal = (-1)^celltype(x,y) * [deriv(x,y,1);deriv(x,y,2)];
                  normal = normal/norm(normal);        % normal n
                  tangent = [-normal(2);normal(1)];    % tangent t
                  kappa = curvature(x,y);          % curvature

                  % mu = mass_dens * nu -> dynamic viscosity
                  mu_2 = visc(celltype(x,y)) * rho_phys(celltype(x,y)); 
                  mu_1 = visc(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))) * rho_phys(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k)));
                  mu_average = (mu_2 + mu_1)*0.5;

                  p_jump = 1/3 * (rho(x,y)*rho_phys(celltype(x,y)) - rho(x+lbm_g.c(1,k),y+lbm_g.c(2,k)) * rho_phys(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))));    
                  % Auf S. 1147 beschrieben, !!! DORT MIT VORZEICHEN-FUCKING-FEHLER
                  % SIEHE doi:10.1016/j.camwa.2009.02.005

                  mu_jump = (mu_2 - mu_1);  % Auch so gewählt, dass außen von innen subtrahiert wird. Ist dann konsistent mit dem Pressure-Jump (mit Vorzeichen-fehler-korrektur! Siehe oben), läuft außerdem stabil (im Gegenteil zu andersrum subtrahiert)

                  % Zwischenergebnisse
                  S_jump_n_n = 1/(2*mu_average) * (p_jump + 2*sigma*kappa) - mu_jump/mu_average * trace(S_average * (normal*normal'));
                  S_jump_n_t = -mu_jump/mu_average * trace(S_average * (normal*tangent')');

                  Lambda_times_S_jump = S_jump_n_n * ((normal'*lbm_g.c(:,k))^2 - (lbm_g.c(:,k)'*lbm_g.c(:,k))/2) + ...
                      2*S_jump_n_t*(normal'*lbm_g.c(:,k))*(tangent'*lbm_g.c(:,k));

                  %% Lambda_i : S^(2)
                  Lambda_times_S_2 = trace(Lambda_i*S_2');

                  %% add_term2 = R_i
                  Lambda_times_A = (-q)*(1-q)*Lambda_times_S_jump - (q-0.5)*Lambda_times_S_2;   % Teilergebnis zur Berechnung von R_i
                  add_term2 = 6*lbm_g.weights(k)*Lambda_times_A;   % R_i

                              % | Don't overwrite interface in opposite cell
                  %% do it    % V Hence, copy cells_new to cells, after inteface-handling done forall cells
                  lbm_g.cells_new(x+lbm_g.c(1,k),y+lbm_g.c(2,k), lbm_g.invDir(k)) = lbm_g.cells(x,y,k) - add_term1 + add_term2;  
                                                                                                     % ^ Minus sign! Draw scenario on paper and note that c_i should be opposite direction (in our interpretation!)
                                                                                                     % | 

              end
        end
    end


    if (max(abs(vel(:))) > 0.1)
      disp('Warning: LBM may become instable, vel > 0.1 [LU]');
    end
    
    data = [data([2:lbm_g.nx-1], [2:lbm_g.ny-1])];

    lbm_g.cells(:,:,:) = lbm_g.cells_new(:,:,:);
    
    % LBM stream
    for i=1:9 
        lbm_g.cells_new(:,:,i) = circshift(lbm_g.cells(:,:,i), [lbm_g.c(1,i),lbm_g.c(2,i),0]); 
    end

    % LBM  ugly (Domain) boundary handling
    for x=1:lbm_g.nx
      %north & South
      lbm_g.cells_new(x  , 2,2) = lbm_g.cells_new(x,1,6);
      lbm_g.cells_new(x  , lbm_g.ny-1,6) = lbm_g.cells_new(x,lbm_g.ny,2) + 6 * lbm_g.weights(2) * lbm_g.c(:,6)'*lbm_g.lidVel;
      if (x<lbm_g.nx)
        lbm_g.cells_new(x+1, 2,3) = lbm_g.cells_new(x,1,7);
        lbm_g.cells_new(x+1, lbm_g.ny-1,5) = lbm_g.cells_new(x,lbm_g.ny,9) + 6 * lbm_g.weights(9) * lbm_g.c(:,5)'*lbm_g.lidVel;
      end
      if (x>1)
        lbm_g.cells_new(x-1, 2,9) = lbm_g.cells_new(x,1,5);
        lbm_g.cells_new(x-1, lbm_g.ny-1,7) = lbm_g.cells_new(x,lbm_g.ny,3) + 6 * lbm_g.weights(3) * lbm_g.c(:,7)'*lbm_g.lidVel;
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

    lbm_g.cells_new(lbm_g.nx, 2:lbm_g.ny-1, :)     = lbm_g.cells_new(2,2:lbm_g.ny-1, :);
    lbm_g.cells_new(1, 2:lbm_g.ny-1, :)  = lbm_g.cells_new(lbm_g.nx-1,2:lbm_g.ny-1, :);

    %copy  cells_new to cells, -> swap old and new buffer
    lbm_g.cells(:,:,:) = lbm_g.cells_new(:,:,:);

    %rho_sum = sum(sum(rho(2:lbm_g.nx-1,2:lbm_g.ny-1)));
    %rho_sum - rho_old   %  >O(1e-10), es wird also irgendwo Masse erzeugt .... :-X ABER: Bekanntes Level-Set Artefakt, siehe Thoemmes Paper
    %rho_old = rho_sum;
  end

  % LBM plot

  figure(4)
  hold off
  plot(sqrt(vel(lbm_g.nx/2,2:lbm_g.ny-1,1).^2+vel(lbm_g.nx/2,2:lbm_g.ny-1,2).^2),[1:lbm_g.ny-2],'ro-');
  hold on

  a2 = 1/(visc(2)/visc(1) + seperator_y*(1-visc(2)/visc(1)));
  a1 = visc(2)/visc(1) * a2;
  offset = a2*seperator_y*(1-visc(2)/visc(1));

  %plot(...
  %  lbm_g.lidVel(1) *...
  %  [a2 * [lbm_g.dx : lbm_g.dx :lbm_g.dx*floor((seperator_y-eps)/lbm_g.dx)],...
  %   (a1 *[lbm_g.dx*ceil(seperator_y/lbm_g.dx) : lbm_g.dx : 1]) + offset],...
  %  [1:lbm_g.ny-2],'b-');
  plot(...
    lbm_g.lidVel(1) *...
    [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*floor((seperator_y-eps)/lbm_g.dx)] -.5*lbm_g.dx),...
    (a1 * ([lbm_g.dx*ceil(seperator_y/lbm_g.dx) : lbm_g.dx : 1]-.5*lbm_g.dx)) + offset],...
    [1:lbm_g.ny-2],'b-');

  plot(([lbm_g.dx:lbm_g.dx:1]-.5*lbm_g.dx)*lbm_g.lidVel(1),[1:lbm_g.ny-2],'g-');
  title('abs(velocity) at x-center column');
  legend({'LBM','analytic','analytic single phase'},'Location','SouthEast')

%  err_rel = norm( (sqrt(vel(lbm_g.nx/2,2:lbm_g.ny-1,1).^2+vel(lbm_g.nx/2,2:lbm_g.ny-1,2).^2)-...
%    (lbm_g.lidVel(1) *...
%    [a2 * [lbm_g.dx : lbm_g.dx :lbm_g.dx*floor((seperator_y-eps)/lbm_g.dx)],...
%     (a1 * [lbm_g.dx*ceil(seperator_y/lbm_g.dx) : lbm_g.dx : 1]) + offset]))...
%     /(lbm_g.lidVel(1) *...
%    [a2 * [lbm_g.dx : lbm_g.dx :lbm_g.dx*floor((seperator_y-eps)/lbm_g.dx)],...
%     (a1 * [lbm_g.dx*ceil(seperator_y/lbm_g.dx) : lbm_g.dx : 1]) + offset]))

  err_rel = norm( (sqrt(vel(lbm_g.nx/2,2:lbm_g.ny-1,1).^2+vel(lbm_g.nx/2,2:lbm_g.ny-1,2).^2)-...
    (lbm_g.lidVel(1) *...
    [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*floor((seperator_y-eps)/lbm_g.dx)] -.5*lbm_g.dx),...
     (a1 * ([lbm_g.dx*ceil(seperator_y/lbm_g.dx) : lbm_g.dx : 1]-.5*lbm_g.dx)) + offset]))...
     /(lbm_g.lidVel(1) *...
    [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*floor((seperator_y-eps)/lbm_g.dx)] -.5*lbm_g.dx),...
     (a1 * ([lbm_g.dx*ceil(seperator_y/lbm_g.dx) : lbm_g.dx : 1]-.5*lbm_g.dx)) + offset]))
     


  err_vek = [err_vek, err_rel];
  figure(5)
  plot(err_vek);
  title(['(Relative) Error in L_2-Norm over \Delta T (:=' num2str(lbm_it) ' LBM-Iterations each)']);

  %Velocity
  figure(2);
  subplot(1,2,1);
  quiver(vel([2:lbm_g.nx-1],[2:lbm_g.ny-1],1)',...
   vel([2:lbm_g.nx-1],[2:lbm_g.ny-1],2)');
  subplot(1,2,2);
  %Density
  %figure(3);
  contourf(rho([2:lbm_g.nx-1],[2:lbm_g.ny-1])')
  title('Density');
  colorbar;


  celltype_old = celltype;


  % "remove" ghost layers
  schemeData.velocity = { lbm_g.dx/lbm_g.dt*vel([2:lbm_g.nx-1], [2:lbm_g.ny-1], 1); lbm_g.dx/lbm_g.dt*vel([2:lbm_g.nx-1], [2:lbm_g.ny-1], 2)};


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
  h = visualizeLevelSet(g, data, 'contour', level, [ 't = ' num2str(tNow) ]);

  % Restore view.
  view(figure_az, figure_el);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);

