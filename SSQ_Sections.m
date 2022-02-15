

temps = [SN1_cell10_fits, SN1_cell11_fits, SN1_cell1_fits, SN1_cell2_fits, SN1_cell3_fits, SN1_cell4_fits, SN1_cell5_fits, SN1_cell6_fits, SN1_cell7_fits, SN1_cell8_fits, SN1_cell9_fits];
for i = 1:length(temps)
COMBO_fits(:,i)=temps(:,i*16); %col.10 is exp 'b'
end

% ctrl_fits = ctrl_fits(:)
% nanmean(nanmean(cell2mat(ctrl_fits)))

%% SSQ
array = {COMBO_cell10, COMBO_cell11, COMBO_cell12, COMBO_cell13, COMBO_cell14, COMBO_cell15, COMBO_cell3, COMBO_cell4, COMBO_cell5, COMBO_cell6, COMBO_cell7, COMBO_cell8};
arrayfits = {COMBO_cell10_fits, COMBO_cell11_fits, COMBO_cell12_fits, COMBO_cell13_fits, COMBO_cell14_fits, COMBO_cell15_fits, COMBO_cell3_fits, COMBO_cell4_fits, COMBO_cell5_fits, COMBO_cell6_fits, COMBO_cell7_fits, COMBO_cell8_fits};

for j = 1:length(array)
    loadcell = array{j};
    loadcellfits = arrayfits{j};
    
for i = 1:16
    if isnan(loadcell{i,8}) == 1
       %ssq(i) = NaN;
       AA1(j,i) = NaN;
    else
%A1 = median(loadcell{i,7}(1:6,2));
%B = loadcell{i,8};
%AA1(i) = abs(A1-B);
AA1(j,i) = abs(loadcellfits{i,16}-loadcell{17,8});
    end
end

%[ssq,~]=sumsqr(AA1);
%deltaN = median(AA1);

end
clear ans B A1 array arrayfits i j loadcell loadcellfits

%%

array ={DOX_SN701, DOX_SN702, DOX_SN703, DOX_SN704, DOX_SN705, DOX_SN706, DOX_SN707, DOX_SN708, DOX_SN709, DOX_SN710};

for i = 1:length(array)
    loadcell = array{i};
    Eapps(i) = loadcell(17,8);
end

clear array loadcell i
