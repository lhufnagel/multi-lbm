function [ type ] = dataToCelltype( data )
%Converts the level set function to celltypes

% Die Funktion geht sehr grob vor, stellt aber sicher, dass es zu keiner
% Interaktion zwischen den beiden Fluiden kommen kann.

env = [1,0;-1,0;0,1;0,-1];

size_data = size(data);

type = ones(size_data);
type_buf = type;

% celltype ohne interface
for i = 1:size_data(1)
    for j = 1:size_data(2)
        if data(i,j) > 0
            type_buf(i,j) = 1;
        elseif data(i,j) < 0
            type_buf(i,j) = -1;
        else
            % data == 0 --> muss interface sein
            type_buf(i,j) = 0;
            continue;
        end
    end
end

type = type_buf;

% % interface
% for i = 2:size_data(1)-1
%     for j = 2:size_data(2)-1
%         %type_sum = sum(sum(type_buf(i-1:i+1,j-1:j+1)));
%         
%         positive = 0;
%         negative = 0;
% 
%         for l = 1:4
%             env_i = type_buf(i+env(l,1),j+env(l,2));
%             if env_i == 1
%                 positive = positive + 1;
%             elseif env_i == -1
%                 negative = negative +1;
%             end
%         end
% 
%         if positive > 0 && negative > 0
%             type(i,j) = 0;
%         else
%             type(i,j) = type_buf(i,j);
%         end
%     end
% 
% end
% 
% % Die Ränder
% for i = 2:size_data(1) - 1
%     % North
%     if (type(1,i+1)+type(1,i-1)) == 0
%         type(1,j) = 0;
%     end
%     % West
%     if (type(i+1,1)+type(i-1,1)) == 0
%         type(j,1) = 0;
%     end
%     % Sout
%     if (type(i+1,size_data(1))+type(i-1,size_data(1))) == 0
%         type(j,size_data(1)) = 0;
%     end
%     % East
%     if (type(size_data(1),i+1)+type(size_data(1),i-1)) == 0
%         type(size_data(1),j) = 0;
%     end
% end