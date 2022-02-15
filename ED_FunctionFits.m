function[] = ED_FunctionFits()
% Function ED_FunctionFits() takes the workspace variables containting 
%   afm_post_pro_editor data and fits exp1 and power2 functions to the
%   post contact E-D raw data and also fits a linear function poly1 
%   to the precontact change in force v depth
%
% Calls user for input: group name, unit name, workspace cell names
% Saves to workspace 3 cell arrays (poly, exp, power)
%   see "Creating Cell Arrays for Data Output" for output array details
%
%
% Robert J. Wiener (c) 2017
%=========================================================================

%afm_post_pro_editor output format
%col.3 = Extension curve (distance-force in 'm' and 'N')
%col.5 = Force curve (depth-force in 'nm' and 'nN')
%col.6 = Pointwise Moduli curve (depth-modulus in 'nm' and 'kPa')


%load workspace cell data
group_name = input('group name ');
unit_name = input('unit name ');
unitgroup_data = input('loadcells ');
unitgroup_data = strsplit(unitgroup_data, ', ');
unitgroup_data = unitgroup_data' ;
for cellnum = 1:length(unitgroup_data)

loadcell = unitgroup_data(cellnum);
loadcell = char(loadcell);
loadcell = evalin('base', loadcell);


%%% INITIALIZE VARIABLES %%%
% % % % % % % % % % % % % % %
%temporary contact point from first cell
%verify pre vector lengths for each cell
nancheck = 1;
while isnan(loadcell{nancheck,8}) == 1
    nancheck = nancheck+1;
end
CP_ind = length(loadcell{nancheck,3}) - length(loadcell{nancheck,5});

%range of precontact depth (-preCD,0)
%range of postcontact depth (0,postCD)
preCD = 800e-9; %e-9 nm
postCD = 200; %pointwise E comes preformatted in 'nm'

% Based on data acq. rate and indentation force
% Number of data points precontact (-800,0 nm) %161
% Number of data points postcontact (0,200 nm) %41
range_pre = CP_ind - length(find(loadcell{nancheck,3}(:,1) < ...
                loadcell{nancheck,3}(CP_ind,1) - preCD))-1;
range_post = length(find(loadcell{nancheck,6}(:,1) < postCD)); %nm
range_post = range_post+1 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%% Begin Data Processing %%%
    for i = 1:length(loadcell)-1
    

%%% Checkpoints %%%
% Checkpoint for deleted data
if isnan(loadcell{i,8}) == 1 
    
% Checkpoint for empty cell    
elseif isempty(loadcell{i,8}) == 1
     precontact_force_data(i,1:range_pre+1) = NaN;
     postcontact_mod_data(i,1:range_post) = NaN;
     
% Checkpoint for not enough precontact data
elseif CP_ind < range_pre
     precontact_force_data(i,1:range_pre+1) = NaN;
     postcontact_mod_data(i,1:range_post) = NaN;
     
else
    
%%% Precontact Data %%%   
%abs distance from contact 'm'
precontact_dist_data(1,:) = -(preCD*1e9):(preCD*1e9)/range_pre:0;
%force 'N'
precontact_force_data(i,:) = loadcell{i,3}(CP_ind - range_pre:CP_ind,2);

%%% Postcontact Data %%%
%distance from contact 'nm'
postcontact_dist_data(1,:) = 0:postCD/range_post:postCD;
%pointwise modulus 'kPa'
postcontact_mod_data(i,:) = loadcell{i,6}(1:range_post,2);
end


% Interpolated pointwise E values at depths
if isempty(loadcell{i,7}) == 0
    InterpE_20 = loadcell{i,7}(1,2);
    InterpE_40 = loadcell{i,7}(3,2);
    InterpE_100 = loadcell{i,7}(9,2);
    InterpE_400 = loadcell{i,7}(39,2);
else
    InterpE_20 = NaN;
    InterpE_40 = NaN;
    InterpE_100 = NaN;
    InterpE_400 = NaN;
end


    end


%Precontact Matrix w/ force-dist data for entire cell
%column1 = abs distance from contact 'm'
%column2 = cell average precontact force
precontact_data(:,1) = precontact_dist_data' ;
for ii = 1:length(precontact_force_data)
    precontact_data(ii,2) = nanmean(precontact_force_data(:,ii));
end
    %NOTE: length(precontact_force_data) == length(loadcell)
    clear ii
    
%Postcontact Matrix w/ E-D data for entire cell
%column1 = distance from contact 'nm'
%column2 = cell average pointwise modulus 'kPa'
postcontact_data(:,1) = postcontact_dist_data' ;
for ii = 1:length(postcontact_mod_data)
    postcontact_data(ii,2) = nanmean(postcontact_mod_data(:,ii));
end
    %NOTE: 
    clear i ii

    
%%% Pre/Post Contact Fucntion Fits %%%
[F_poly1,poly_error] = fit(precontact_data(:,1),precontact_data(:,2),'poly1');
[F_exp1,exp_error] = fit(postcontact_data(2:end,1),postcontact_data(2:end,2),'exp1');
[F_power2,power_error] = fit(postcontact_data(2:end,1),postcontact_data(2:end,2),'power2');



%%% Creating Cell Arrays for Data Output %%%

%Precontact cell array
% (1) Group
% (2) Unit
% (3) Precontact Data
% (4) Polynomial Fit Coefficient 'p1'
% (5) Polynomial Fit Coefficient 'p2'  %f(x)= p1x + p2
% (6) Rsquared
precontactFits{cellnum,1} = group_name;
precontactFits{cellnum,2} = unit_name;
precontactFits{cellnum,3} = precontact_data;
precontactFits{cellnum,4} = F_poly1.p1;
precontactFits{cellnum,5} = F_poly1.p2;
precontactFits{cellnum,6} = poly_error.rsquare;

%Postcontact cell array exp1
% (1) Group
% (2) Unit
% (3) Postcontact Data
% (4) average pointwise cell E @ 20nm
% (5) average pointwise cell E @ 40nm
% (6) average pointwise cell E @ 100nm
% (7) average pointwise cell E @ 400nm
% (8) exp1 Fit Coefficient 'a'
% (9) exp1 Fit Coefficient 'b'  %f(x)= a*exp(b*x)
% (10) Rsquared
postcontactFits_exp1{cellnum,1} = group_name;
postcontactFits_exp1{cellnum,2} = unit_name;
postcontactFits_exp1{cellnum,3} = postcontact_data;
postcontactFits_exp1{cellnum,4} = InterpE_20;
postcontactFits_exp1{cellnum,5} = InterpE_40;
postcontactFits_exp1{cellnum,6} = InterpE_100;
postcontactFits_exp1{cellnum,7} = InterpE_400;
postcontactFits_exp1{cellnum,8} = F_exp1.a;
postcontactFits_exp1{cellnum,9} = F_exp1.b;
postcontactFits_exp1{cellnum,10} = exp_error.rsquare;

%Postcontact cell array power2
% (1) Group
% (2) Unit
% (3) Postcontact Data
% (4) average pointwise cell E @ 20nm
% (5) average pointwise cell E @ 40nm
% (6) average pointwise cell E @ 100nm
% (7) average pointwise cell E @ 400nm
% (8) power2 Fit Coefficient 'a'
% (9) power2 Fit Coefficient 'b'
% (10) power2 Fit Coefficient 'c'  %f(x)= a*x^b+c
% (11) Rsquared
postcontactFits_power2{cellnum,1} = group_name;
postcontactFits_power2{cellnum,2} = unit_name;
postcontactFits_power2{cellnum,3} = postcontact_data;
postcontactFits_power2{cellnum,4} = InterpE_20;
postcontactFits_power2{cellnum,5} = InterpE_40;
postcontactFits_power2{cellnum,6} = InterpE_100;
postcontactFits_power2{cellnum,7} = InterpE_400;
postcontactFits_power2{cellnum,8} = F_power2.a;
postcontactFits_power2{cellnum,9} = F_power2.b;
postcontactFits_power2{cellnum,10} = F_power2.c;
postcontactFits_power2{cellnum,11} = power_error.rsquare;


% clear variables
clear loadcell cellnum loadcell CP_ind
clear postcontact_dist_data postcontact_mod_data
clear precontact_dist_data precontact_force_data
clear postcontact_data postcontact_dist_data postcontact_mod_data 
clear precontact_force_data range_post range_pre precontact_data precontact_dist_data
clear CP_ind exp_error F_exp1 F_poly1 F_power2 poly_error power_error
clear InterpE_100 InterpE_20 InterpE_40 InterpE_400  
clear postCD preCD 

end


%save cell arrays with cell names
assignin('base', horzcat(unit_name, 'poly1'), precontactFits);
assignin('base', horzcat(unit_name, 'exp1'), postcontactFits_exp1);
assignin('base', horzcat(unit_name, 'power2'), postcontactFits_power2);





   
% % % % -------------------------------------------------------- % % %    
% %%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot_opt = 0; % % % % % % % %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
% if plot_opt == 1
% % plot raw extension curve from contact point
% figure
% plot(loadcell{1,3}(contact_point_ind:end,1),loadcell{1,3}(contact_point_ind:end,2))


end

