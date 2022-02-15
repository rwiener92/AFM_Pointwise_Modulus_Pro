
ctrl_cells = ctrl_800;
%endoBq_cells = ENDBq_100;
%endo_cells = ENDO_100;
%combo_cells = COMBO_100;
SN1_cells = SN1_800;
SN2_cells = SN2_800;
SN2Bq_cells = SN2Bq_800;

[ctrl_H,ctrl_B] = hist(ctrl_cells,20);
ctrl_HN = ctrl_H ./ sum(ctrl_H);
% [enbq_H,enbq_B] = hist(endoBq_cells,ctrl_B);
% enbq_HN = enbq_H ./ sum(enbq_H);
% [endo_H,endo_B] = hist(endo_cells,ctrl_B);
% endo_HN = endo_H ./ sum(endo_H);
% [combo_H,combo_B] = hist(combo_cells,ctrl_B);
% combo_HN = combo_H ./ sum(combo_H);

[SN1_H,SN1_B] = hist(SN1_cells,ctrl_B);
SN1_HN = SN1_H ./ sum(SN1_H);
[SN2_H,SN2_B] = hist(SN2_cells,ctrl_B);
SN2_HN = SN2_H ./ sum(SN2_H);
[SN2Bq_H,SN2Bq_B] = hist(SN2Bq_cells,ctrl_B);
SN2Bq_HN = SN2Bq_H ./ sum(SN2Bq_H);

figure
% plot(ctrl_B,ctrl_HN,'-k',ctrl_B,combo_HN,'-r',ctrl_B,endo_HN,'-b',ctrl_B,enbq_HN,'-g',...
%         ctrl_B,SN1_HN,'y',ctrl_B,SN2_HN,'m',ctrl_B,SN2Bq_HN,'c');
%legend({'CTRL';'HYAL + Hep';'ENDO';'ENDO + Bq';'SN1';'SN2';'SN2 + Bq'});
plot(ctrl_B,ctrl_HN,'-k', ctrl_B,SN1_HN,'y',ctrl_B,SN2_HN,'m',ctrl_B,SN2Bq_HN,'c');
legend({'CTRL';'SN1';'SN2';'SN2 + Bq'});
xlabel('Brush Layer Resistance  ( N/m )')
ylabel('Normalized Frequency')
title('Precontact (-800,0 nm)')

clear
