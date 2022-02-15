
NIL_cell1{37,8} = nanmean([NIL_cell1{1:36,8}]);
NIL_cell2{37,8} = nanmean([NIL_cell2{1:36,8}]);
NIL_cell3{37,8} = nanmean([NIL_cell3{1:36,8}]);
NIL_cell4{37,8} = nanmean([NIL_cell4{1:36,8}]);
NIL_cell5{37,8} = nanmean([NIL_cell5{1:36,8}]);
NIL_cell6{37,8} = nanmean([NIL_cell6{1:36,8}]);
NIL_cell7{37,8} = nanmean([NIL_cell7{1:36,8}]);
NIL_cell8{37,8} = nanmean([NIL_cell8{1:36,8}]);
NIL_cell9{37,8} = nanmean([NIL_cell9{1:36,8}]);
NIL_cell10{37,8} = nanmean([NIL_cell10{1:36,8}]);
NIL_cell11{37,8} = nanmean([NIL_cell11{1:36,8}]);
NIL_cell12{37,8} = nanmean([NIL_cell12{1:36,8}]);
NIL_cell13{37,8} = nanmean([NIL_cell13{1:36,8}]);
NIL_cell14{37,8} = nanmean([NIL_cell14{1:36,8}]);



NIL_fit_avgs=[NIL_cell1(37,8), NIL_cell2(37,8), NIL_cell3(37,8), NIL_cell4(37,8), NIL_cell5(37,8), NIL_cell6(37,8), NIL_cell7(37,8), NIL_cell8(37,8), NIL_cell9(37,8), NIL_cell10(37,8), NIL_cell11(37,8), NIL_cell12(37,8), NIL_cell13(37,8), NIL_cell14(37,8)];

NIL_avg_ptw=[NIL_cell1{37,7}(:), NIL_cell2{37,7}(:), NIL_cell3{37,7}(:), NIL_cell4{37,7}(:), NIL_cell5{37,7}(:), NIL_cell6{37,7}(:), NIL_cell7{37,7}(:), NIL_cell8{37,7}(:), NIL_cell9{37,7}(:), NIL_cell10{37,7}(:), NIL_cell11{37,7}(:), NIL_cell12{37,7}(:), NIL_cell13{37,7}(:), NIL_cell14{37,7}(:)];

NIL_fit_avgs = cell2mat(NIL_fit_avgs);
NIL_fit_avgs(imag(NIL_fit_avgs) ~= 0) = NaN;

NIL_avg_ptw(imag(NIL_avg_ptw) ~= 0) = NaN;

ptw_NIL(:,1) = mean(NIL_avg_ptw,2);
ptw_NIL(:,2) = std(NIL_avg_ptw,0,2)/sqrt(length(NIL_fit_avgs'));


figure
errorbar(depth(1:49),ptw_ctrl(1:49,1),ptw_ctrl(1:49,2),'ko')
hold on; plot([0,depth(49)],[mean(ctrl_fit_avgs), mean(ctrl_fit_avgs)],'k')

hold on

% errorbar(depth(1:49),ptw_IMA(1:49,1),ptw_IMA(1:49,2),'go')
% hold on; plot([0,depth(49)],[mean(IMA_fit_avgs), mean(IMA_fit_avgs)],'g')
% 
% hold on

errorbar(depth(1:49),ptw_VAN(1:49,1),ptw_VAN(1:49,2),'bo')
hold on; plot([0,depth(49)],[mean(VAN_fit_avgs), mean(VAN_fit_avgs)],'b')

hold on

errorbar(depth(1:49),ptw_DAS(1:49,1),ptw_DAS(1:49,2),'ro')
hold on; plot([0,depth(49)],[mean(DAS_fit_avgs), mean(DAS_fit_avgs)],'r')


xlabel('Indentation Depth, Z (nm)')
ylabel('Elastic Modulus (kPa)')
legend( 'Control Pointwise','Control Apparent',...
    'Vandatinib Pointwise','Vandatinib Apparent',...
    'Dasatinib Pointwise','Dasatinib Apparent' )


