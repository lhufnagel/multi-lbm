%Lid driven cavity on [0,x_len] x [0, y_len] for t in [0, t_end], top lid has velocity lidVel in x-Direction
%D2Q9 LBM, BGK, Stream-Collide. SiWiR-2 Inspired ;)

t_end = 5; % [s]
x_len = 1; % [m]
y_len = 1; % [m]
lidVel = 2; % [m/s]
visc  = 1e-4;% [m^2/s]

g=Grid;
g.dx=0.02; % [m]
g.dt=0.001; % [s]
g.c_s=sqrt(1/3); % [m/s]
g.lidVel = [lidVel*g.dt/g.dx; 0];
g.omega = 1/(visc/(g.c_s^2*g.dt) + 1/2);
g.nx = x_len/g.dx + 2; %Ghost layer
g.ny = y_len/g.dx + 2; %Ghost layer
g.type=celltype.Regular;

%init
g.cells = ones(g.nx,g.ny,9);
for i=1:9
  g.cells(:,:,i)=g.weights(i);
end

g.type([1:g.nx], [1:g.ny]) = celltype.Regular;

%go for it
for t=1:(t_end/g.dt)

  %stream
  for i=1:9 
      g.cells_new(:,:,i) = circshift(g.cells(:,:,i), [g.c(1,i),g.c(2,i),0]); 
  end

  %ugly boundary handling
  for x=1:g.nx
    %north & South
    g.cells_new(x  , 2,2) = g.cells_new(x,1,6);
    g.cells_new(x  , g.ny-1,6) = g.cells_new(x,g.ny,2) - 2/(g.c_s^2) * g.weights(2) * g.c(:,2)'*g.lidVel;
    if (x<g.nx)
      g.cells_new(x+1, 2,3) = g.cells_new(x,1,7);
      g.cells_new(x+1, g.ny-1,5) = g.cells_new(x,g.ny,9) - 2/(g.c_s^2) * g.weights(9) * g.c(:,9)'*g.lidVel;
    end
    if (x>1)
      g.cells_new(x-1, 2,9) = g.cells_new(x,1,5);
      g.cells_new(x-1, g.ny-1,7) = g.cells_new(x,g.ny,3) - 2/(g.c_s^2) * g.weights(3) * g.c(:,3)'*g.lidVel;
    end
  end

  for y=2:g.ny-1
    %east & west

    %solid boundaries
   g.cells_new(g.nx-1, y ,8) = g.cells_new(g.nx,y, 4);
   g.cells_new(g.nx-1, y-1, 7) = g.cells_new(g.nx,y, 3);
   g.cells_new(g.nx-1, y+1, 9) = g.cells_new(g.nx,y, 5);
   g.cells_new(2, y  ,4) = g.cells_new(1,y, 8);
   g.cells_new(2, y-1,5) = g.cells_new(1,y, 9);
   g.cells_new(2, y+1,3) = g.cells_new(1,y, 7);

    %periodic boundaries
    %g.cells_new(2, y, 4)     = g.cells_new(g.nx,y, 4);
    %g.cells_new(g.nx-1, y,8) = g.cells_new(1,y,8);
    %if (y>2)
    %  g.cells_new(2, y, 3)     = g.cells_new(g.nx,y, 3);
    %  g.cells_new(g.nx-1, y,9) = g.cells_new(1,y,9);
    %end
    %if (y<g.ny-1)
    %  g.cells_new(2, y,5)     = g.cells_new(g.nx,y, 5);
    %  g.cells_new(g.nx-1, y,7) = g.cells_new(1,y,7);
    %end
  end


  %collide
  rho = zeros(g.nx,g.ny);
  vel = zeros(g.nx,g.ny,2);

  for i=1:9
    rho(:,:) = rho(:,:) + g.cells_new(:,:,i);
    vel(:,:,1) = vel(:,:,1) + g.c(1,i) * g.cells_new(:,:,i);
    vel(:,:,2) = vel(:,:,2) + g.c(2,i) * g.cells_new(:,:,i);
  end

  %! Here cells_new is implicitly swapped with the old cells -> Stream-Collide!
  for i=1:9
    cTimesU = g.c(1,i) * vel(:,:,1) + g.c(2,i) * vel(:,:,2);
    g.cells(:,:,i) = g.cells_new(:,:,i) - g.omega .* (g.cells_new(:,:,i) - ...
        g.weights(i) .* (rho(:,:) + ...
          1/(g.c_s^2) .* cTimesU(:,:) + ...
          1/(2*g.c_s^4) .* (cTimesU(:,:)).^2 - ...
          1/(2*g.c_s^2) .* (vel(:,:,1).^2 + vel(:,:,2).^2)));
  end

  fprintf('Steps: %d \tDone: %d\n',(t_end/g.dt),t)

  %plot
  if (mod(t,25)==0)
    %fig = figure(1);
    %set(fig,'Position', [0 0 800 300]);

    figure(1);
    subplot(1,2,1);
    quiver(vel([2:g.nx-1],[2:g.ny-1],1)',...
     vel([2:g.nx-1],[2:g.ny-1],2)');
    subplot(1,2,2);
    contourf(rho([2:g.nx-1],[2:g.ny-1])')
    colorbar;
  end
end
