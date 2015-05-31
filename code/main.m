%D2Q9 LBM
%Lid driven cavity on [0,x_len] x [0, y_len] for t in [0, t_end]

t_end = 5; % [s]
x_len = 1; % [m]
y_len = 2; % [m]
visc  = 1e-6;% [m^2/s]

g=grid;
g.dx=0.1;
g.dt=0.1;
g.c_s=sqrt(1/3);
g.omega = 1/(visc/(g.c_s*g.c_s*g.dt) + 1/2);
g.nx = x_len/g.dt;
g.ny = y_len/g.dt;
g.cells = zell.empty(g.nx, 0); 

%init
for x=1:g.nx
  for y=1:g.ny
    g.cells(x, y) = zell;

    g.cells(x,y).type = celltype.Regular;
    g.cells(x,y).pdfs = g.weights;

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
      g.cells(x,y).type = celltype.NorthSolid;
    end

  end
end

for t=1:(t_end/g.dt)

  for x=1:g.nx
    for y=1:g.ny
      g.cells(x, y).stream(g, x, y);
    end
  end

  for x=1:g.nx
    for y=1:g.ny
      g.cells(x, y).collide(g);
    end
  end

end
