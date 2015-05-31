classdef Grid
  properties
    dx
    dt
    nx
    ny
    omega
    c_s
    lidVel
    cells

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
  end
end
