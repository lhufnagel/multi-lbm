clear;
% Remark: Origin (1,1) is bottom left!
%
%   ^
% y |__>
%     x

t_end = 10; % [s]
x_len = 1; % [m]
y_len = 1; % [m] %Ignoriert von Level-Set Bibilothek!
rho_phys(1) = 1000; % [kg/m^3] % Mass-density of Fluid 1, e.g. ~1000 for Water, 1,2 for Air. In Lattice-Boltzmann-Units Cell-density is varying around 1!
rho_phys(2) = 1; % [kg/m^3] %
lidVel = 0; % [m/s]
visc(1) = 16;% [m^2/s] kinematic viscosity of fluid 1
visc(2) = 16;% [m^2/s] kinematic viscosity of fluid 2
sigma = 1.6e-8;  % [kg/s^2] surface tension the two phases, e.g. ~ 76*10^-3 between water and air
lbm_it = 10; % Number of LBM-iterations until level-set update

% Diese Zahl sollte mMn folgendermaﬂen beschr‰nkt sein:
% Delta_T = lbm_it * lbm_g.dt ist das Zeitinterval, in dem sich Level-Set und LBM abwechseln.
% Die maximale Geschwindigkeit in der Lid-Driven-Cavity ist lidVel (wenn man starke Kr¸mmungseffekte und daraus folgende groﬂe Oberfl‰chenspannungen an kleinen Blasen ignoriert).
% Da ich mˆchte, dass sich das Interface MAXIMAL eine Zelle weit bewegt pro Delta_T (CFL-Bedingung, auﬂerdem geht sonst der Refill-Algorithmus kaputt),
% muss lidVel * Delta_T = lidVel * lbm_it * lbm_g.dt < lbm_g.dx sein!

lbm_g=LBM_Grid;
lbm_g.dx=0.1*x_len; % [m]
lbm_g.dt=lbm_g.dx^2;% [s]  %Diffusive scaling..


lbm_g.lidVel = lbm_g.dt/lbm_g.dx * [lidVel; 0];
lbm_g.omega(1) = 1/(3*visc(1)*lbm_g.dt/lbm_g.dx^2 + 1/2);
lbm_g.omega(2) = 1/(3*visc(2)*lbm_g.dt/lbm_g.dx^2 + 1/2);

if (max(lbm_g.omega) > 1.7) %|| min(lbm_g.omega) < 0.5)  %%TODO abkl‰ren? bound von unten: arXiv:0812.3242v2 letzter Absatz 1/1.8 = 0.555
  disp('WARNING: omega out of range; interface handling maybe instable! Change dx, dt or viscosity');
end

lbm_g.nx = round(x_len/lbm_g.dx) + 2; %Ghost layer
lbm_g.ny = round(y_len/lbm_g.dx) + 2; %Ghost layer

%init LBM data
lbm_g.cells = ones(lbm_g.nx,lbm_g.ny,9);
for i=1:9
  lbm_g.cells(:,:,i) = lbm_g.weights(i);
end

%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
run('./addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
plotSteps = t_end/(lbm_g.dt*lbm_it) + 1;         % How many intermediate plots to produce?

% Period at which intermediate plots should be produced.
tPlot = t_end / (plotSteps - 1);

%---------------------------------------------------------------------------
% Pause after each plot?
pauseAfterPlot = 0;
% Delete previous Interface plot before showing next?
deleteLastPlot = 1;

%---------------------------------------------------------------------------
% Use periodic boundary conditions?
level_set_periodic = 1;

% Create the Level-Set-grid.
lsm_g.dim = 2;
lsm_g.min = 0.5*lbm_g.dx;
lsm_g.dx = lbm_g.dx;
if(level_set_periodic)
  lsm_g.max = (x_len - 0.5*lsm_g.dx);
  lsm_g.bdry = @addGhostPeriodic;
else
  lsm_g.max = (x_len - 0.5*lsm_g.dx);
  lsm_g.bdry = @addGhostExtrapolate;
end
lsm_g = processGrid(lsm_g);

%---------------------------------------------------------------------------
% Create initial conditions (a circle/sphere).
%   Note that in the level_set_periodic BC case, these initial conditions will not
%   be continuous across the boundary unless the circle is perfectly centered.
%   In practice, we'll just ignore that little detail.

%   data         Implicit surface function at t_max.
%   Lsm_g        Level-Set-Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

testcase = 'circle';

switch(testcase)
  case 'circle'
    center = [ 0.5; .5]*x_len;
    radius = 0.3*x_len;
    
    data = zeros(size(lsm_g.xs{1}));
    data = data + (lsm_g.xs{1} - center(1)).^2 + (lsm_g.xs{2} - center(2)).^2;
    data = sqrt(data) - radius;
  case 'line'
    % Horizontale Trennlinie
    seperator_y = 0.42*x_len;
    data = lsm_g.xs{2} - seperator_y;
  otherwise
    error('Testcase does not exist: %s', testcase);
end
data0 = data;

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

accuracy = 'veryHigh';
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

%---------------------------------------------------------------------------
% Initialize Level Set Display
f = figure(1);
hold off;

h = visualizeLevelSet(lsm_g, data, 'contour', 0, [ 't = ' num2str(0) ]);

hold on;
axis(lsm_g.axis);
daspect([ 1 1 1 ]);

%---------------------------------------------------------------------------
% Loop until t_end (subject to a little roundoff).
tNow = 0;

celltype_old = dataToCelltype(data);
celltype_old = [celltype_old(:,1), celltype_old, celltype_old(:,end)];
celltype_old = [celltype_old(end,:); celltype_old; celltype_old(1,:)];

err_vek = [];
mass_vek = [];
pressure_vek =[];

startTime = cputime;

while(t_end - tNow > 100 * eps * t_end)
  
  % Pass Data from level set to LBM
  
  %Contour line of the zero level set
  C=contourc(lsm_g.xs{1}(:,1), lsm_g.xs{2}(1,:), data, [0 0]);
  C=C(:,2:end);
  
  % get curvature
  [curvature, ~] = curvatureSecond(lsm_g, data);
  
  % get first order derivative (=> normals)
  deriv = zeros(lbm_g.nx-2,lbm_g.ny-2,2);
  [derivL,derivR] = upwindFirstENO3(lsm_g,data,1);
  deriv(:,:,1) = 0.5 * (derivL + derivR);
  [derivL,derivR] = upwindFirstENO3(lsm_g,data,2);
  deriv(:,:,2) = 0.5 * (derivL + derivR);
  
  % "Adding" Ghost layers here!
  % Easier accessing with LBM indices
  data = [data(:,1), data, data(:,end)];
  data = [data(end,:); data; data(1,:)];
  
  deriv = [deriv(:,1,:), deriv, deriv(:,end,:)];
  deriv = [deriv(end,:,:); deriv; deriv(1,:,:)];
  
  curvature = [curvature(:,1), curvature, curvature(:,end)];
  curvature = [curvature(end,:); curvature; curvature(1,:)];
  
  % Refill-Algorithm for cells that changed their type
  celltype = dataToCelltype(data);
  changed_type = celltype - celltype_old;
  [x, y] = find(changed_type);
  
  for i=1:length(x)
    if x(i) < 3 || x(i) > lsm_g.N(1) || y(i) < 3 || y(i) > lsm_g.N(2)
      continue;
    end
    
    normal = [deriv(x(i),y(i),1);deriv(x(i),y(i),2)];
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
    
    % rho and vel should be up-to-date from previous LBM (before Level-Set)
    
    rho_interp = 3 * rho(x(i) + c_imax(1), y(i) + c_imax(2)) - ...
      3 * rho(x(i) + 2*c_imax(1), y(i) + 2*c_imax(2)) + ...
      rho(x(i) + 3*c_imax(1), y(i) + 3*c_imax(2));
    
    % TODO Evtl aus n‰chstem Timestep interpolieren?!
    % UND: Ich verstehe nicht, wieso die im Paper den vektor abziehen... er zeigt doch nach innen und ich interpoliere von inneren Punkten..
    
    vel_interp = 3 * vel(x(i) + c_imax(1), y(i) + c_imax(2),:) - ...
      3 * vel(x(i) + 2*c_imax(1), y(i) + 2*c_imax(2),:) + ...
      vel(x(i) + 3*c_imax(1), y(i) + 3*c_imax(2),:);
    
    % TODO hier wirds sehr unschˆn: In der Interpolationsformel im Thˆmmes-Paper wird f¸r u_Gamma die Geschwindigkeit v aus dem Level-Set (?!?! Siehe auch Paper Kapitel 5.3) verwendet.
    % Diese Geschwindigkeit stammt - so wie ich es verstehe - aus der Fast-Marching/Constant extension velocity  methode.
    % Diese Methode implementiert die Bibliothek von Mitchel allerdings nicht. (Im Gegensatz zu evtl http://ktchu.serendipityresearch.org/software/lsmlib/ )
    % -> ??
    % ... -> ich verwende die Formel von rho auch f¸r die Geschwindigkeit.
    % Die Formel ist auch zu finden in Lallemand & Luo: Lattice Boltzmann method for moving boundaries. J Comput Phys 184(2): 406- 421
    % und sollte laut dem Paper ‰hnliche Ergebnisse produzieren (siehe Ende von Kapitel 3 in dem Paper)
    
    
    equil=zeros(1,9);
    non_eq=zeros(1,9);

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
    
    %sync cells_new to cells, as we use it temporarily here
    lbm_g.cells_new(:,:,:) = lbm_g.cells(:,:,:);
    
    %Interface behandlung PRE-stream
    
    for x=1:lbm_g.nx
      for y = 1:lbm_g.ny
        for k = 2:9
          
          if x+lbm_g.c(1,k) < 1 || x+lbm_g.c(1,k) > lbm_g.nx || y+lbm_g.c(2,k) < 1 || y+lbm_g.c(2,k) > lbm_g.ny
            continue
          end
          
          % Hier evtl die "isNearInterface(..)"-Funktion verwenden? (Doku S. 124)
          % Das ist mit die teuerste Abfrage..
          if celltype(x,y) == celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))
            continue
          end
          
          
          % q
          q = data(x+lbm_g.c(1,k), y+lbm_g.c(2,k))/(data(x+lbm_g.c(1,k), y+lbm_g.c(2,k)) - data(x,y));
          
          % Die Liniensegmente vom "exakten" Inteface kriegt man mit contourc(lsm_g.xs{1}(:,1), lsm_g.xs{2}(1,:), data, [0 0]);
          % R¸ckgabe-Format Doku: http://de.mathworks.com/help/matlab/ref/contour-properties.html#prop_ContourMatrix
          % Eventuell kann man damit q noch genauer berechnen (also exakt den Schnittpunkt mit dem Link c_i bilden)
          % N.B.: Mit dem contourc-Aufruf arbeitet zumindest der visualizeLevelSet-Aus der Mitchel-Bibliothek intern
          
          if (x > 2 && x < lbm_g.nx-1 && y > 2 && y < lbm_g.ny-1)
            x_2 = [lsm_g.xs{1}(x-1,y-1);lsm_g.xs{2}(x-1,y-1)];%die stimmen hier noch nicht!
            x_1 = [lsm_g.xs{1}(x+lbm_g.c(1,k)-1,y+lbm_g.c(2,k)-1); lsm_g.xs{2}(x+lbm_g.c(1,k)-1,y+lbm_g.c(2,k)-1)];
            
            [x_int,y_int]=intersections([x_1(1),x_2(1)],[x_1(2),x_2(2)],C(2,:),C(1,:),false);
            if (~isempty(x_int))
              q = norm(x_1-[x_int(1);y_int(1)])/norm(x_1-x_2);
            end
          end
          
          
          
          
          
          
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
          % Normale, Tangente und Kr¸mmung werden in der Toolbox bestimmt
          normal = (-1)^celltype(x,y) * [deriv(x,y,1);deriv(x,y,2)];
          normal = normal/norm(normal);        % normal n
          tangent = [-normal(2);normal(1)];    % tangent t
          kappa = curvature(x,y);          % curvature
          
          % mu = mass_dens * nu -> dynamic viscosity
          mu_2 = visc(celltype(x,y)) * rho_phys(celltype(x,y));
          mu_1 = visc(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))) * rho_phys(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k)));
          mu_average = (mu_2 + mu_1)*0.5;
          
          p_jump = 1/3 *(rho(x,y)*rho_phys(celltype(x,y)) - rho(x+lbm_g.c(1,k),y+lbm_g.c(2,k)) * rho_phys(celltype(x+lbm_g.c(1,k),y+lbm_g.c(2,k))));
          % Auf S. 1147 beschrieben, !!! DORT MIT VORZEICHEN-FUCKING-FEHLER
          % SIEHE doi:10.1016/j.camwa.2009.02.005
          
          mu_jump = (mu_2 - mu_1);  % Auch so gew‰hlt, dass auﬂen von innen subtrahiert wird. Ist dann konsistent mit dem Pressure-Jump (mit Vorzeichen-fehler-korrektur! Siehe oben), l‰uft auﬂerdem stabil (im Gegenteil zu andersrum subtrahiert)
          
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
    
    %Copy back the (PRE-Stream) interface-treated links from cells_new to cells
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
      
      %%periodic boundaries
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
    
    %Guarantee periodic boundaryies,
    lbm_g.cells_new(lbm_g.nx, 2:lbm_g.ny-1, :)     = lbm_g.cells_new(2,2:lbm_g.ny-1, :);
    lbm_g.cells_new(1, 2:lbm_g.ny-1, :)  = lbm_g.cells_new(lbm_g.nx-1,2:lbm_g.ny-1, :);
    
    %copy  cells_new to cells, -> swap old and new buffer
    lbm_g.cells(:,:,:) = lbm_g.cells_new(:,:,:);
    
    
    p1 = 1/3*rho_phys(1)*mean(rho(celltype==1));
    p2 = 1/3*rho_phys(2)*mean(rho(celltype==2));
    pressure_vek = [pressure_vek, p1-p2];
    
    %rho_sum = sum(sum(rho(2:lbm_g.nx-1,2:lbm_g.ny-1)));
    %mass_vek = [mass_vek, rho_sum];
  end
  
  rho = zeros(lbm_g.nx,lbm_g.ny);
  vel = zeros(lbm_g.nx,lbm_g.ny,2);
  
  for i=1:9
    rho(:,:) = rho(:,:) + lbm_g.cells(:,:,i);
    vel(:,:,1) = vel(:,:,1) + lbm_g.c(1,i) * lbm_g.cells(:,:,i);
    vel(:,:,2) = vel(:,:,2) + lbm_g.c(2,i) * lbm_g.cells(:,:,i);
  end
  
  % LBM plot
  
  %Velocity
  figure(2);
  subplot(1,2,1);
  %Scale vectors, such that, 0.05 velocity in [LU] := 1 length
  quiver(1/0.05*vel( 2:lbm_g.nx-1 , 2:lbm_g.ny-1 ,1)',...
    1/0.05*vel(2:lbm_g.nx-1 ,2:lbm_g.ny-1,2)','AutoScale','off');
  
  title(['Velocity, Range: [0, ' num2str(max(vel(:))) ']']);
  subplot(1,2,2);
  %Density
  contourf(rho(2:lbm_g.nx-1,2:lbm_g.ny-1)')
  title('Density');
  colorbar;
  
  % Mass losses
  % figure(3);
  % plot(mass_vek/((lbm_g.nx-2)*(lbm_g.ny-2))-ones(size(mass_vek)));
  % title(['\Delta Mass (relative) over iterations']);
  
  % Error
  % figure(4)
  % hold off
  % plot(sqrt(vel(lbm_g.nx/2,2:lbm_g.ny-1,1).^2+vel(lbm_g.nx/2,2:lbm_g.ny-1,2).^2),[1:lbm_g.ny-2],'ro-');
  % hold on
  %
  % a2 = (1+lbm_g.dx)/(visc(2)*rho_phys(2)/(visc(1)*rho_phys(1)) + seperator_y*(1-rho_phys(2)*visc(2)/(visc(1)*rho_phys(1))));
  % a1 = rho_phys(2)*visc(2)/(visc(1)*rho_phys(1)) * a2;
  % offset = a2*seperator_y*(1-visc(2)*rho_phys(2)/(visc(1)*rho_phys(1)))-lbm_g.dx;
  %
  % plot(...
  %   lbm_g.lidVel(1) *...
  %   [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*ceil((seperator_y-eps)/lbm_g.dx)] -.5*lbm_g.dx),...
  %   (a1 * ([lbm_g.dx*ceil(1+seperator_y/lbm_g.dx-eps) : lbm_g.dx : 1]-.5*lbm_g.dx)) + offset],...
  %   [1:lbm_g.ny-2],'b-');
  
  % err_rel = norm( abs(sqrt(vel(lbm_g.nx/2,2:lbm_g.ny-1,1).^2+vel(lbm_g.nx/2,2:lbm_g.ny-1,2).^2)-...
  %   (lbm_g.lidVel(1) *...
  %   [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*ceil((seperator_y-eps)/lbm_g.dx)] -.5*lbm_g.dx),...
  %   (a1 * ([lbm_g.dx*ceil(1+seperator_y/lbm_g.dx-eps) : lbm_g.dx : 1]-.5*lbm_g.dx)) + offset]))...
  %    /(lbm_g.lidVel(1) *...
  %   [a2 * ([lbm_g.dx : lbm_g.dx :lbm_g.dx*ceil((seperator_y-eps)/lbm_g.dx)] -.5*lbm_g.dx),...
  %   (a1 * ([lbm_g.dx*ceil(1+seperator_y/lbm_g.dx-eps) : lbm_g.dx : 1]-.5*lbm_g.dx)) + offset]))
  
  % err_vek = [err_vek, err_rel];
  % xlabel('Velocity [LU]');
  % ylabel('y, Cell index');
  % title('abs(velocity) at x-center column');
  % legend({'LBM','analytic'},'Location','SouthEast')
  
  % figure(5)
  % plot(err_vek);
  % title(['(Relative) Error in L_2-Norm over \Delta T (:=' num2str(lbm_it) ' LBM-Iterations each)']);
  
  figure(5)
  plot(pressure_vek);
  hold on
  pressure_jump_analytic = 2*sigma/radius;
  plot(pressure_jump_analytic*ones(size(pressure_vek)),'g:');
  hold off
  title('Pressure jump between interior and exterior of bubble');
  %axis([0, length(pressure_vek), 0,2*pressure_jump_analytic]);
  
  
  celltype_old = celltype;
  % "remove" ghost layers
  schemeData.velocity = { lbm_g.dx/lbm_g.dt*vel(2:lbm_g.nx-1, 2:lbm_g.ny-1, 1); lbm_g.dx/lbm_g.dt*vel(2:lbm_g.nx-1, 2:lbm_g.ny-1, 2)};
  data = data(2:lbm_g.nx-1, 2:lbm_g.ny-1);
  
  %% Level Set code
  
  % Reshape data array into column vector for ode solver call.
  y0 = data(:);
  
  % How far to step?
  tSpan = [ tNow, min(t_end, tNow + tPlot) ];
  
  % Take a timestep.
  [ t, y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
    integratorOptions, schemeData);
  tNow = t(end);
  
  % Get back the correctly shaped data array
  data = reshape(y, lsm_g.shape);
  
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
  
  % Create new visualization.
  h = visualizeLevelSet(lsm_g, data, 'contour', 0, [ 't = ' num2str(tNow) ]);
  
  % Restore view.
  view(figure_az, figure_el);
  
end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);

