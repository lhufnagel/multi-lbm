classdef LBM_Grid
  properties
    dx
    dt
    nx
    ny
    omega % Use 2-dim-Vector.... a little ugly
    lidVel % die gibt es nur einmal
    cells
    cells_new

    % D2Q9; Indexing:
    % 9   2   3
    %     ^
    % 8 < 1 > 4
    %     v
    % 7   6   5
    %
    % Remark: Origin (1,1) is bottom left!
    %
    %   ^
    % y |__>
    %     x

    weights = [ 4/9,   1/9,   1/36,  1/9,   1/36,   1/9,    1/36,    1/9,    1/36];
    c       = [[0;0], [0;1], [1;1], [1;0], [1;-1], [0;-1], [-1;-1], [-1;0], [-1;1]];
    c_len   = sqrt(sum([[0;0], [0;1], [1;1], [1;0], [1;-1], [0;-1], [-1;-1], [-1;0], [-1;1]].^2));
    invDir  = [0, 6, 7, 8, 9, 2, 3, 4, 5];

  end
end
