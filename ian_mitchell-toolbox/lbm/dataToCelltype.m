function [ type ] = dataToCelltype( data )
%Converts the level set function to celltypes

type = sign(data);
type(type == 0) = 1;
type = 1.5 - 0.5*sign(type);

% type contains 1's where data is >= 0
% and 2's where data < 0. 
% Hence, normals generated from the gradient of data will always point towards the 1's!

% "Adding" Ghost layers here! Limited to 2D..
type = [type(:,1), type, type(:,end)];
type = [type(end,:); type; type(1,:)];
