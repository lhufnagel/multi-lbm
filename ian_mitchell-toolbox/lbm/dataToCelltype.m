function [ type ] = dataToCelltype( data )
%Converts the level set function to celltypes

type = sign(data);
type(type == 0) = 1;
type = 1.5 - 0.5*sign(type);

% "Adding" Ghost layers here! Limited to 2D..
type = [type(:,end), type, type(:,1)];
type = [type(1,:); type; type(end,:)];
