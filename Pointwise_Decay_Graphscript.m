

val1 = {ctrl_cell9, ctrl_cell8, ctrl_cell6, ctrl_cell5, ctrl_cell4, ctrl_cell3, ctrl_cell2, ctrl_cell14, ctrl_cell13, ctrl_cell12, ctrl_cell11, ctrl_cell10, ctrl_cell1};

val2 = {COMBO_cell9, COMBO_cell8, COMBO_cell7, COMBO_cell6, COMBO_cell5, COMBO_cell4, COMBO_cell3, COMBO_cell2, COMBO_cell15, COMBO_cell14, COMBO_cell13, COMBO_cell12, COMBO_cell11, COMBO_cell10};

val3 = {};
val4 = {};
val5 = {};
val6 = {};

% full depth of pointwise indentation data (x-axis)
X = 20:10:400;

figure

if isempty(val1) == 0
[~,a(1,1)] = size(val1);
for i_1 = 1:a(1)
   
   plot(X, val1{1,i_1}{17,7},'r')
   hold on
            
end
end

if isempty(val2) == 0
[~,a(1,2)] = size(val2);
for i_2 = 1:a(2)
   
   hold on
   plot(X, val2{1,i_2}{17,7},'b')
            
end
end

if isempty(val3) == 0
[~,a(1,3)] = size(val3);
for i_3 = 1:a(3)
   
   hold on
   plot(X, val3{1,i_3}{17,7},'g')
            
end
end

if isempty(val4) == 0
[~,a(1,4)] = size(val4);
for i_4 = 1:a(4)
   
   hold on
   plot(X, val4{1,i_4}{17,7},'m')
            
end
end

if isempty(val5) == 0
[~,a(1,5)] = size(val5);
for i_5 = 1:a(5)
   
   hold on
   plot(X, val5{1,i_5}{17,7},'y')
            
end
end

if isempty(val6) == 0
[~,a(1,6)] = size(val6);
for i_6 = 1:a(6)
   
   hold on
   plot(X, val6{1,i_6}{17,7},'c')
            
end
end
    
clear a X
clear i_1 i_2 i_3 i_4 i_5 i_6
clear val1 val2 val3 val4 val5 val6

title('Part2')
legend('r-ctrl','b-COMBO')
       