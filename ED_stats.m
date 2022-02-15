
cell_load = {part9_161111exp1, part9_161111poly1, part9_161111power2};
cell1 = cell_load{2}; %poly
cell2 = cell_load{1}; %exp
cell3 = cell_load{3}; %power

%poly cell{2}
p1=[mean(cell2mat(cell1(:,4))),median(cell2mat(cell1(:,4))),var(cell2mat(cell1(:,4)))];
p2=[mean(cell2mat(cell1(:,5))),median(cell2mat(cell1(:,5))),var(cell2mat(cell1(:,5)))];
r2_poly=mean(cell2mat(cell1(:,6)));

%exp cell{1}
a_exp=[mean(cell2mat(cell2(:,8))),median(cell2mat(cell2(:,8))),var(cell2mat(cell2(:,8)))];
b_exp=[mean(cell2mat(cell2(:,9))),median(cell2mat(cell2(:,9))),var(cell2mat(cell2(:,9)))];
r2_exp=mean(cell2mat(cell2(:,10)));
exp_depth = [nanmean(cell2mat(cell2(:,4))),nanmean(cell2mat(cell2(:,5))),nanmean(cell2mat(cell2(:,6))),nanmean(cell2mat(cell2(:,7)))];

%power cell{3}
a_power=[mean(cell2mat(cell3(:,8))),median(cell2mat(cell3(:,8))),var(cell2mat(cell3(:,8)))];
b_power=[mean(cell2mat(cell3(:,9))),median(cell2mat(cell3(:,9))),var(cell2mat(cell3(:,9)))];
c_power=[mean(cell2mat(cell3(:,10))),median(cell2mat(cell3(:,10))),var(cell2mat(cell3(:,10)))];
r2_power=mean(cell2mat(cell3(:,11)));
power_depth = [nanmean(cell2mat(cell3(:,4))),nanmean(cell2mat(cell3(:,5))),nanmean(cell2mat(cell3(:,6))),nanmean(cell2mat(cell3(:,7)))];
