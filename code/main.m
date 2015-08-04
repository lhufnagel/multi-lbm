%Lid driven cavity on [0,x_len] x [0, y_len] for t in [0, t_end], top lid has velocity lidVel in x-Direction
%D2Q9 LBM, BGK, Stream-Collide. SiWiR-2 Inspired ;)

t_end = 5; % [s]
x_len = 1; % [m]
y_len = 1; % [m]
lidVel = 0.3; % [m/s]
visc  = 1e-3;% [m^2/s]

g=Grid;
g.dx=0.05; % [m]
g.dt=0.01; % [s]
g.c_s=sqrt(1/3); % [m/s]
g.lidVel = [lidVel*g.dt/g.dx; 0];
g.omega = 1/(visc/(g.c_s^2*g.dt) + 1/2);
g.nx = x_len/g.dx + 2; %Ghost layer
g.ny = y_len/g.dx + 2; %Ghost layer
g.cells = Cell;

%init
for x=1:g.nx
  for y=1:g.ny
    g.cells(x, y) = Cell;

    g.cells(x,y).type = celltype.Regular;
    g.cells(x,y).pdfs = g.weights;

    %Below ordering is critical! Do better not touch...
    if (x == 1)
      g.cells(x,y).type = celltype.WestPeriodic;
    end
    if (x == g.nx)
      g.cells(x,y).type = celltype.EastPeriodic;    
    end
    if (y == 1)
      g.cells(x,y).type = celltype.SouthSolid;
    end
    if (y == g.ny)
      g.cells(x,y).type = celltype.NorthMovSolid;
    end
    
  end
end

densField = ones(g.nx-2,g.ny-2);
velField_u = zeros(g.nx-2,g.ny-2);
velField_v = zeros(g.nx-2,g.ny-2);

fig = figure(1);
set(fig,'Position', [0 0 800 300]);

%go for it
for t=1:(t_end/g.dt)
  for x=1:g.nx
    for y=1:g.ny
      g.cells(x, y).stream(g, x, y);
    end
  end

  for y=1:g.ny
      g.cells(1, y).boundaryHandling(g, 1, y);
      g.cells(g.nx, y).boundaryHandling(g, g.nx, y);
  end
  for x=1:g.nx
      g.cells(x, 1).boundaryHandling(g, x, 1);
      g.cells(x, g.ny).boundaryHandling(g, x, g.ny);
  end


  for x=1:g.nx
    for y=1:g.ny
      [rho, v] = g.cells(x, y).collide(g);
      if (x > 1 && x <g.nx &&...
        y > 1 && y < g.ny)
        densField(x-1,y-1)=rho;
        velField_u(x-1,y-1)=v(1);
        velField_v(x-1,y-1)=v(2);
      end
    end
  end

  fprintf('Steps: %d \tDone: %d\n',(t_end/g.dt),t)
  if (mod(t,2)==0)
    figure(1);
    subplot(1,2,1);
    quiver(velField_u',velField_v');
    subplot(1,2,2);
    contourf(densField')
  end
end
