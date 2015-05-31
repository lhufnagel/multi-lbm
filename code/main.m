%Lid driven cavity on [0,x_len] x [0, y_len] for t in [0, t_end], top lid has velocity lidVel in x-Direction
%D2Q9 LBM, BGK, Stream-Collide. SiWiR-2 Insipred ;)

t_end = 5; % [s]
x_len = 1; % [m]
y_len = 1; % [m]
lidVel = 0.2; % [m/s]
visc  = 1e-2;% [m^2/s]

g=Grid;
g.dx=0.05; % [m]
g.dt=0.01; % [s]
g.c_s=sqrt(1/3); % [m/s]
g.lidVel = [lidVel*g.dt/g.dx; 0];
g.omega = 1/(visc/(g.c_s^2*g.dt) + 1/2);
g.nx = x_len/g.dx;
g.ny = y_len/g.dx;
g.cells = Cell;

%init
for x=1:g.nx
  for y=1:g.ny
    g.cells(x, y) = Cell;

    g.cells(x,y).type = celltype.Regular;
    g.cells(x,y).pdfs = g.weights;
    g.cells(x,y).pdfs_new = g.weights;

    if (x == 1)
      g.cells(x,y).type = celltype.WestSolid;
    end

    if (x == g.nx)
      g.cells(x,y).type = celltype.EastSolid;
    end

    if (y == 1)
      g.cells(x,y).type = celltype.SouthSolid;
    end

    if (y == g.ny)
      g.cells(x,y).type = celltype.NorthMovSolid;
    end

  end
end

densField = ones(g.nx,g.ny);
velField_u = zeros(g.nx,g.ny);
velField_v = zeros(g.nx,g.ny);

%go for it
for t=1:(t_end/g.dt)
  for x=1:g.nx
    for y=1:g.ny
      g.cells(x, y).stream(g, x, y);
    end
  end

  for x=1:g.nx
    for y=1:g.ny
      g.cells(x, y).boundaryHandling(g, x, y);
    end
  end


  for x=1:g.nx
    for y=1:g.ny
      [rho, v] = g.cells(x, y).collide(g);
      densField(x,y)=rho;
      velField_u(x,y)=v(1);
      velField_v(x,y)=v(2);
    end
  end

  disp(t)
  if (mod(t,2)==0)
    figure(1);
    quiver(velField_u',velField_v');
    figure(2);
    contourf(densField')
  end
end
