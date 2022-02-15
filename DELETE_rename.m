% Script changes 'interpolated' cell into cell of NaN
% SECTION B: Recalculates average interpolatied depth dependent modulus in final row 
cellname = ctrl_cell1;

% For all 16 force curves interpolated data that was 'Deleted' become NaN
check_opt = 0;
for i = 1:36
%     if isequaln(cellname{i,8},NaN) == 1
%        temp = cellname{i,7};
%        temp(:,:) = NaN; 
%        cellname{i,7} = temp;
%        check_opt = 1;

    if cellname{i,8}>25
    cellname{i,8} = NaN(79,2);
    cellname{i,7} = NaN(79,2);
    check_opt = 1;
    end
    
end

            %if check_opt == 1
%%%%% SECTION B %%%%%   
 for ii = 1:36
    E_temp = cellname{ii,7};
    mean_temp(:,ii) = E_temp(:,2);
 end
 
 for iii = 1:39
    final_mean(iii,1) = nanmean(mean_temp(iii,:));
 end
    cellname{37,7} = final_mean;
    
    clear i ii iii
    clear mean_temp temp final_mean E_temp
        
ctrl_cell1 = cellname;

clear cellname check_opt
            %else
            %end

 