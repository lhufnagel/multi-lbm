function [ type ] = dataToCelltype( data )
%Converts the level set function to celltypes

size_data = size(data);

type = ones(size_data);

% celltype ohne interface
for i = 1:size_data(1)
    for j = 1:size_data(2)
        if data(i,j) < 0
            type(i,j) = -1;
        end
    end
end
