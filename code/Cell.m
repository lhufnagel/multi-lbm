classdef Cell < handle
  properties
    % D2Q9; Indexing:
    % 9   2   3
    %     ^
    % 8 < 1 > 4
    %     v
    % 7   6   5
    %
    %   ^
    % y |__>
    %     x

    pdfs    = [ 0,      0,     0,     0,      0,      0,      0,      0,      0];
    pdfs_new= [ 0,      0,     0,     0,      0,      0,      0,      0,      0];
    type;
  end
  methods
    function stream(obj, g, x, y)
      if (obj.type ~= celltype.Regular)
        return
      end

      for i=1:9
        g.cells(x + g.c(1,i), y + g.c(2,i)).pdfs_new(i) = obj.pdfs(i);
      end
    end

    function boundaryHandling(obj, g, x, y)
      switch (obj.type)
        case celltype.EastSolid
          g.cells(x-1, y  ).pdfs_new(8) = obj.pdfs_new(4);
          g.cells(x-1, y-1).pdfs_new(7) = obj.pdfs_new(3);
          g.cells(x-1, y+1).pdfs_new(9) = obj.pdfs_new(5);
        case celltype.WestSolid
          g.cells(x+1, y  ).pdfs_new(4) = obj.pdfs_new(8);
          g.cells(x+1, y-1).pdfs_new(5) = obj.pdfs_new(9);
          g.cells(x+1, y+1).pdfs_new(3) = obj.pdfs_new(7);
        case celltype.SouthSolid
          g.cells(x  , y+1).pdfs_new(2) = obj.pdfs_new(6);
          if (x<g.nx)
            g.cells(x+1, y+1).pdfs_new(3) = obj.pdfs_new(7);
          end
          if (x>1)
            g.cells(x-1, y+1).pdfs_new(9) = obj.pdfs_new(5);
          end
        case celltype.NorthMovSolid
          g.cells(x  , y-1).pdfs_new(6) = obj.pdfs_new(2) - 2/(g.c_s^2) * g.weights(2) * g.c(:,2)'*g.lidVel;
          if (x<g.nx)
            g.cells(x+1, y-1).pdfs_new(5) = obj.pdfs_new(9) - 2/(g.c_s^2) * g.weights(9) * g.c(:,9)'*g.lidVel;
          end
          if (x>1)
            g.cells(x-1, y-1).pdfs_new(7) = obj.pdfs_new(3) - 2/(g.c_s^2) * g.weights(3) * g.c(:,3)'*g.lidVel;
          end
      end
    end

    function [rho, vel] = collide(obj, g)
      %Macroscopic measures [LU]
      rho = 1.;
      vel = [0; 0];

      if (obj.type ~= celltype.Regular)
        return;
      end

      rho=0.;

      for i=1:9
        rho = rho + obj.pdfs_new(i);
        vel = vel + obj.pdfs_new(i) * g.c(:,i);
      end

      %! Here the array pdfs_new is implicitly swapped with the old pdfs -> Stream-Collide!
      for i=1:9
        obj.pdfs(i) = obj.pdfs_new(i) - g.omega * (obj.pdfs_new(i) - g.weights(i)*(rho + ...
          1/(g.c_s^2) * g.c(:,i)' * vel + ...
          1/(2*g.c_s^4) * (g.c(:,i)' * vel)^2 - ...
          1/(2*g.c_s^2) * (vel' * vel)^2));
      end
    end
  end
end
