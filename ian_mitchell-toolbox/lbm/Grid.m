classdef Grid
  properties
    dx
    dt
    nx
    ny
    rho_1
    rho_2
    omega_1
    omega_2
    c_s
    lidVel_1
    lidVel_2
    cells
    cells_new
    %type

    % D2Q9; Indexing:
    % 9   2   3
    %     ^
    % 8 < 1 > 4
    %     v
    % 7   6   5
    %
    % 
    % Remark: Origin (1,1) is bottom left!
    %
    %   ^
    % y |__>
    %     x

    weights = [ 4/9,   1/9,   1/36,  1/9,   1/36,   1/9,    1/36,    1/9,    1/36];
    c       = [[0;0], [0;1], [1;1], [1;0], [1;-1], [0;-1], [-1;-1], [-1;0], [-1;1]];
    invDir  = [0, 6, 7, 8, 9, 2, 3, 4, 5];

  end
  methods
    function boundaryHandling(obj, g, x, y) %not working, pass by value issue.........:-X
      switch (g.type(x,y))
        case celltype.EastSolid
          g.cells_new(x-1, y ,8) = obj.cells_new(x,y, 4);
          g.cells_new(x-1, y-1, 7) = obj.cells_new(x,y, 3);
          g.cells_new(x-1, y+1, 9) = obj.cells_new(x,y, 5);
        case celltype.EastPeriodic
          g.cells_new(2, y, 4)     = obj.cells_new(x,y, 4);
          if (y>2)
            g.cells_new(2, y, 3)     = obj.cells_new(x,y, 3);
          end
          if (y<g.ny-1)
            g.cells_new(2, y,5)     = obj.cells_new(x,y, 5);
          end
        case celltype.WestSolid
          g.cells_new(x+1, y  ,4) = obj.cells_new(x,y, 8);
          g.cells_new(x+1, y-1,5) = obj.cells_new(x,y, 9);
          g.cells_new(x+1, y+1,3) = obj.cells_new(x,y, 7);
        case celltype.WestPeriodic
          g.cells_new(g.nx-1, y,8) = obj.cells_new(x,y,8);
          if (y>2)
            g.cells_new(g.nx-1, y,9) = obj.cells_new(x,y,9);
          end
          if (y<g.ny-1)
            g.cells_new(g.nx-1, y,7) = obj.cells_new(x,y,7);
          end
        case celltype.SouthSolid
          g.cells_new(x  , y+1,2) = obj.cells_new(x,y,6);
          if (x<g.nx)
            g.cells_new(x+1, y+1,3) = obj.cells_new(x,y,7);
          end
          if (x>1)
            g.cells_new(x-1, y+1,9) = obj.cells_new(x,y,5);
          end
        case celltype.NorthMovSolid
          g.cells_new(x  , y-1,6) = obj.cells_new(x,y,2) - 2/(g.c_s^2) * g.weights(2) * g.c(:,2)'*g.lidVel;
          if (x<g.nx)
            g.cells_new(x+1, y-1,5) = obj.cells_new(x,y,9) - 2/(g.c_s^2) * g.weights(9) * g.c(:,9)'*g.lidVel;
          end
          if (x>1)
            g.cells_new(x-1, y-1,7) = obj.cells_new(x,y,3) - 2/(g.c_s^2) * g.weights(3) * g.c(:,3)'*g.lidVel;
          end
      end
    end
  end
end
