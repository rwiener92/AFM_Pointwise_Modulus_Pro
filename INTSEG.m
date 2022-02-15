%INTSEG creates boxplots for all cell types at specified depth intervals
% Part9
% ###_INTSEG is cell array
% X-axis is cell#
% Y-axis is depth point [20 40 60 80 100, 320 340 360 380 400]
% Column 7 is a vertical contact of that entire row

DI = [1:2:9, 31:2:39];
for ii = 1:10
    
for iii = 1:16
    
ctrl_cell_INTSEG{ii,1}(iii,1) = ctrl_cell1{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,2}(iii,1) = ctrl_cell2{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,3}(iii,1) = ctrl_cell3{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,4}(iii,1) = ctrl_cell4{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,5}(iii,1) = ctrl_cell5{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,6}(iii,1) = ctrl_cell6{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,7}(iii,1) = ctrl_cell7{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,8}(iii,1) = ctrl_cell8{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,9}(iii,1) = ctrl_cell9{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,10}(iii,1) = ctrl_cell10{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,11}(iii,1) = ctrl_cell11{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,12}(iii,1) = ctrl_cell12{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,13}(iii,1) = ctrl_cell13{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,14}(iii,1) = ctrl_cell14{iii,7}(DI(ii),2);
ctrl_cell_INTSEG{ii,15}(iii,1) = ctrl_cell15{iii,7}(DI(ii),2);

A_cell_INTSEG{ii,1}(iii,1) = ENDO_cell1{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,2}(iii,1) = ENDO_cell2{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,3}(iii,1) = ENDO_cell3{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,4}(iii,1) = ENDO_cell4{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,5}(iii,1) = ENDO_cell5{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,6}(iii,1) = ENDO_cell7{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,7}(iii,1) = ENDO_cell8{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,8}(iii,1) = ENDO_cell9{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,9}(iii,1) = ENDO_cell10{iii,7}(DI(ii),2);
A_cell_INTSEG{ii,10}(iii,1) = ENDO_cell11{iii,7}(DI(ii),2);

B_cell_INTSEG{ii,1}(iii,1) = COMBO_cell1{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,2}(iii,1) = COMBO_cell3{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,3}(iii,1) = COMBO_cell4{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,4}(iii,1) = COMBO_cell5{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,5}(iii,1) = COMBO_cell6{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,6}(iii,1) = COMBO_cell7{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,7}(iii,1) = COMBO_cell8{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,8}(iii,1) = COMBO_cell10{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,9}(iii,1) = COMBO_cell11{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,10}(iii,1) = COMBO_cell12{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,11}(iii,1) = COMBO_cell13{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,12}(iii,1) = COMBO_cell14{iii,7}(DI(ii),2);
B_cell_INTSEG{ii,14}(iii,1) = COMBO_cell15{iii,7}(DI(ii),2);

C_cell_INTSEG{ii,1}(iii,1) = SN2_cell1{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,2}(iii,1) = SN2_cell2{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,3}(iii,1) = SN2_cell3{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,4}(iii,1) = SN2_cell4{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,5}(iii,1) = SN2_cell5{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,6}(iii,1) = SN2_cell6{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,7}(iii,1) = SN2_cell7{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,8}(iii,1) = SN2_cell8{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,9}(iii,1) = SN2_cell9{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,10}(iii,1) = SN2_cell10{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,11}(iii,1) = SN2_cell11{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,12}(iii,1) = SN2_cell12{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,13}(iii,1) = SN2_cell13{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,14}(iii,1) = SN2_cell14{iii,7}(DI(ii),2);
C_cell_INTSEG{ii,15}(iii,1) = SN2_cell15{iii,7}(DI(ii),2);

D_cell_INTSEG{ii,1}(iii,1) = SN1_cell1{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,2}(iii,1) = SN1_cell2{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,3}(iii,1) = SN1_cell3{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,4}(iii,1) = SN1_cell5{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,5}(iii,1) = SN1_cell7{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,6}(iii,1) = SN1_cell8{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,7}(iii,1) = SN1_cell9{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,8}(iii,1) = SN1_cell11{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,9}(iii,1) = SN1_cell12{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,10}(iii,1) = SN1_cell13{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,11}(iii,1) = SN1_cell14{iii,7}(DI(ii),2);
D_cell_INTSEG{ii,12}(iii,1) = SN1_cell15{iii,7}(DI(ii),2);
end

ctrl_cell_INTSEG{ii,7} = cat(1, ctrl_cell_INTSEG{ii,1:6});
A_cell_INTSEG{ii,7} = cat(1, A_cell_INTSEG{ii,1:6});
B_cell_INTSEG{ii,7} = cat(1, B_cell_INTSEG{ii,1:6});
C_cell_INTSEG{ii,7} = cat(1, C_cell_INTSEG{ii,1:6});
D_cell_INTSEG{ii,7} = cat(1, D_cell_INTSEG{ii,1:6});
end

%Box Plots
figure
% 2x5 subplots

depth_title = [20 40 60 80 100 320 340 360 380 400];
for ij = 1:10
    
subplot(2,5,ij)
boxplot([ctrl_cell_INTSEG{ij,7},A_cell_INTSEG{ij,7},B_cell_INTSEG{ij,7}...
    C_cell_INTSEG{ij,7},D_cell_INTSEG{ij,7}],...
            'Labels',{'ctrl','Endo','COMBO','SN2','SN1'})
ylim([0 100])
title(horzcat(num2str(depth_title(ij)),' nm'))  
end

clear ii iii ij
clear DI depth_title
clear ctrl_cell_INTSEG A_cell_INTSEG B_cell_INTSEG C_cell_INTSEG D_cell_INTSEG




