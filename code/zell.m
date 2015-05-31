% cell as name is already taken -.-

classdef zell
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
    type;
  end
  methods
    function stream(obj, g, x, y)
      if (obj.type == celltype.Regular)
        for i=1:9
          g.cells(x + g.c(1,i), y + g.c(2,i)).pdfs(i) = obj.pdfs(i);
        end
      end
    end

    function [rho, vel] = collide(obj, g);
      %Macroscopic measures [LU]
      rho = 0.;
      vel = [0; 0];

      for i=1:9
        rho = rho + obj.pdfs(i);
        vel = vel + obj.pdfs(i) * g.c(:,i)
      end

    end
  end
end
