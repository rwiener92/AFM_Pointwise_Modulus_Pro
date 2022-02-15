%function[] = ED_FunctionFits_full()
% Function ED_FunctionFits() takes the workspace variables containing
%   afm_post_pro_editor data and fits exp1 and power2 functions to the
%   post contact E-D raw data and also fits a linear function poly1 
%   to the precontact change in force v depth
%
% NOTE: Unlike the ED_FunctionFits, the full version returns a single
%        "cell_fits" array with all FM data points per cell
%
% Calls user for input: group name, unit name, workspace cell names
% Saves to workspace cell array containting: (poly, exp, power)
%   see "Creating Cell Arrays for Data Output" for output array details
%
%
% Robert J. Wiener (c) 2017
%=========================================================================

%afm_post_pro_editor output format
%col.3 = Extension curve (distance-force in 'm' and 'N')
%col.5 = Force curve (depth-force in 'nm' and 'nN')
%col.6 = Pointwise Moduli curve (depth-modulus in 'nm' and 'kPa')

%%== RUN AS SCRIPT ==%%
clear ans
vars = struct2cell(whos)' ;
[vars_runs,~]=size(vars);

for i = 1:vars_runs

%%== RUN AS FN ==%%    
%load workspace cell data
% group_name = input('group name ');
% unit_name = input('unit name ');
% unitgroup_data = input('loadcells ');
group_name = vars(i);
unit_name = horzcat(vars(i), char(i));
unitgroup_data = char(vars(i));
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
medEDind = find(loadcell{3,7}(:,1)==50); %50nm depth

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
    precontact_force_data(i,1:range_pre+1) = NaN;
    postcontact_mod_data(i,1:range_post) = NaN;
    
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
    
        clear i
        for i = 1:length(loadcell)-1

%disp(unitgroup_data(cellnum))   
%Precontact Matrix w/ force-dist data for entire Force Map
%column1 = abs distance from contact 'm'
%column2 = cell precontact force
precontact_data(:,1) = precontact_dist_data(1,:)' ;
precontact_data(:,2) = precontact_force_data(i,:)' ;
    %NOTE: length(precontact_force_data) == length(loadcell)
    
%Postcontact Matrix w/ E-D data for entire cell
%column1 = distance from contact 'nm'
%column2 = cell average pointwise modulus 'kPa'
postcontact_data(:,1) = postcontact_dist_data(1,3:end)' ;
postcontact_data(:,2) = postcontact_mod_data(i,2:end)' ;


c = isnan(precontact_data(:,2));
c = sum(c);
if c > 0
    disp(unitgroup_data(cellnum)) 
else
%%% Pre/Post Contact Fucntion Fits %%%
[F_poly1,poly_error] = fit(precontact_data(:,1),precontact_data(:,2),'poly1');
[F_exp1,exp_error] = fit(postcontact_data(:,1),postcontact_data(:,2),'exp1');
[F_power2,power_error] = fit(postcontact_data(:,1),postcontact_data(:,2),'power2');


%%% Creating Cell Arrays for Data Output %%%

% Cell Array Fit Data
% (1) Group
% (2) X-coordinate of Indentation
% (3) Y-coordinate of Indentation
% (4) Precontact Data
% (5) Postcontact Data
% (6) Polynomial Fit Coefficient 'p1' 
% (7) Polynomial Fit Coefficient 'p2'  %f(x)= p1x + p2
% (8) poly Rsquared
% (9) exp1 Fit Coefficient 'a'
% (10) exp1 Fit Coefficient 'b'  %f(x)= a*exp(b*x)
% (11) exp Rsquared
% (12) power2 Fit Coefficient 'a'
% (13) power2 Fit Coefficient 'b'
% (14) power2 Fit Coefficient 'c'  %f(x)= a*x^b+c
% (15) power Rsquared
% (16) Median ED %depth range = medEDind
CellFits{i,1} = horzcat(group_name, unit_name);
CellFits{i,2} = cell2mat(loadcell(i,1));
CellFits{i,3} = cell2mat(loadcell(i,2));
CellFits{i,4} = precontact_data ;
CellFits{i,5} = postcontact_data ;
CellFits{i,6} = F_poly1.p1;
CellFits{i,7} = F_poly1.p2;
CellFits{i,8} = poly_error.rsquare;
CellFits{i,9} = F_exp1.a;
CellFits{i,10} = F_exp1.b;
CellFits{i,11} = exp_error.rsquare;
CellFits{i,12} = F_power2.a;
CellFits{i,13} = F_power2.b;
CellFits{i,14} = F_power2.c;
CellFits{i,15} = power_error.rsquare;
CellFits{i,16} = median(loadcell{3,7}(1:medEDind,2));


end
%%% Checkpoints %%%
% fit checkpoint: kill numbers with bad Rsquared
if c>0
  CellFits(i,1:15) = {NaN};
end

if CellFits{i,8} < 0.5
	CellFits{i,6} = NaN ;
    CellFits{i,7} = NaN ;

elseif CellFits{i,11} < 0.7
	CellFits{i,9} = NaN ;
    CellFits{i,10} = NaN ;

elseif CellFits{i,15} < 0.7
	CellFits{i,12} = NaN ;
    CellFits{i,13} = NaN ;
    CellFits{i,14} = NaN ;

end


        end
           
        
%save cell arrays with cell names
assignin('base', horzcat(char(unitgroup_data(cellnum)), '_fits'), CellFits);

% clear variables
clear loadcell cellnum CP_ind
clear postcontact_dist_data postcontact_mod_data
clear precontact_dist_data precontact_force_data
clear postcontact_data postcontact_dist_data postcontact_mod_data 
clear precontact_force_data range_post range_pre precontact_data precontact_dist_data
clear CP_ind exp_error F_exp1 F_poly1 F_power2 poly_error power_error
clear InterpE_100 InterpE_20 InterpE_40 InterpE_400  
clear postCD preCD 
end


   
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

disp('Program Completed')
clear i unit_name unitgroup_data vars vars_runs group_name
clear c CellFits medEDind nancheck


