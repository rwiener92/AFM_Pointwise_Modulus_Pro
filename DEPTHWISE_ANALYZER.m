function data_summary = DEPTHWISE_ANALYZER(indenter,C_S)
% MODULUS_SUMMARY = DEPTHWISE_ANALYZER(IND_TYPE,STIFFNESS) returns the
% pointwise elastic moduli from 300 to 2000 nm indentation depths with 50nm
% increments for a given force curve with a FV001LQ file. IND_TYPE is an
% integer designator for the indenter geometry where the indenter could be
% one of the following:
%       1- Standard blunt cone
%       2- Custom sphere
%       3- Custom hemisphere
%
% STIFFNESS is the cantilever spring coefficient in N/m.
% 
% AZELOGLU E.U. (C) 2008

global min_depth max_depth stp_depth I_R
% Minimum and maximum depths and the step increments required for
% interpolation of pointwise modulus:
min_depth = 100;
max_depth = 1200;
stp_depth = 100; % in nanometers
% Indenter radius (if a hemispherical or spherical tip is being used:
I_R = 12.5e-6; % in meters
% Load a vector of indentation depths to be interpolated for:
depth_vector = (min_depth:stp_depth:max_depth)';

% Ask for a debug list feed which has the list of indentation numbers on
% the first column and a numerical indicator of 0 (automatically accept
% file) or 1 (ask the user to debug) on the second column:
[file_temp,path_name] =  uigetfile({'*.csv','CSV Files (*.csv)',...
        ; '*.*',  'All Files (*.*)'}, ...
        'Select the feed with pre-processed data list');
if file_temp == 0
    clear file_temp
    clear path_name
    clear sigma
    return
end
debug_temp = csvread([path_name '/' file_temp]);
ind_temp = debug_temp(:,1); log_temp = debug_temp(:,2);
debug_indents = ind_temp(logical(log_temp));
clear *_temp


% Ask for a directory of result files to be processed
current_directory = pwd;
path_name = uigetdir(current_directory,'Choose Directory');
path_name = [path_name '/'];
% Load the summary file(s) listing into a directory listing vector
file_list = dir([path_name '*.fv001lq']);
% Kill program if the directory listing returns an empty handle
if isempty(file_list) == 1
    warning('Requested file(s) cannot be found!');
    data_summary = [];
    return
end


% DATA GROUPING
% Ask the user for the data groups and the beginning and end indentation
% numbers that envelope the force curve array:
prompt = {'First Array','First Indent','Last Indent',...
    'Second Array','First Indent','Last Indent',...
    'Third Array','First Indent','Last Indent',...
    'Fourt Array','First Indent','Last Indent'};
dlg_title = 'FC Array Information';
num_lines= 1;
def     = {'samp1','1','25',...
    'samp2','26','50',...
    'samp3','51','75',...
    'samp4','76','100'};
options.Resize='on';
options.WindowStyle='normal';
answer  = inputdlg(prompt,dlg_title,num_lines,def,options);
if isempty(answer) == 1
    clear dlg_title num_lines def answer prompt
    return
end     
% Create a 2-by-4 matrix that has the first and last indentation for each
% data group that's being sorted. First column is the beginning and end of
% the first array
indent_range(1,1) = str2double(answer{02});
indent_range(2,1) = str2double(answer{03});
indent_range(1,2) = str2double(answer{05});
indent_range(2,2) = str2double(answer{06});
indent_range(1,3) = str2double(answer{08});
indent_range(2,3) = str2double(answer{09});
indent_range(1,4) = str2double(answer{11});
indent_range(2,4) = str2double(answer{12});
clear dlg_title num_lines def prompt


% Recognize the indentation numbers from the file name and convert them
% into integers from text:
for k = 1:length(file_list)
    name_str = file_list(k).name; n = 1;
    while n <= length(name_str)
        if strcmp(name_str(1,n),'.') == 1
            indent_temp = name_str(n+1:n+3);
            indent_temp = str2double(indent_temp);
            break
        end
        n = n + 1;
    end
    indent_no(k) = indent_temp;
    clear indent_temp
end

% ------------------------
% ANALYSIS LOOP
% Go over each FC array (four groups) and load the indentation data into a
% vector, which will then be either automatically analyzed (the curve will
% be considered good and the contact point will be accepted as is) or fed
% to user for debugging (either change contact point or reject curve all
% together.
%
% Note: Last modulus number in the results vector is a "fitted modulus"
% value obtained by least-squares fitting the force-depth vector. Currently
% this is only done for the spherical and hemispherical indentations.
for k = 1:4
   data_summary{1,k} = answer{3*k-2};
   data_summary{2,k} = depth_vector;
   m = 1;
   for n = 1:length(file_list)
       if indent_no(n) >= indent_range(1,k) && ...
               indent_no(n) <= indent_range(2,k)
           if isempty(find(debug_indents == indent_no(n), 1)) == 1
               modulus(:,m) = auto_accept_curve(...
                   path_name,indent_no(n),indenter,C_S);
               edit_list(1,m) = indent_no(n);
               edit_list(2,m) = 0;
               disp(['Completed: ' data_summary{1,k} ' curve ' ...
                   num2str(indent_no(n)) ' of ' num2str(indent_no(end))]);
               m = m + 1;
           else
               modulus(:,m) = user_accept_curve(...
                   path_name,indent_no(n),indenter,C_S);
               edit_list(1,m) = indent_no(n);
               edit_list(2,m) = 1;
               disp(['Completed: ' data_summary{1,k} ' curve ' ...
                   num2str(indent_no(n)) ' of ' num2str(indent_no(end))]);
               m = m + 1;
           end
       end
   end
   data_summary{3,k} = modulus;
   data_summary{4,k} = edit_list;
   clear modulus
end
% The final result is presented as a cell array as follows:
% Row 1 = Name of the grouped samples
% Row 2 = Depth vector for pointwise analysis used for the given group
% Row 3 = Pointwise modulus as a function of depth vectors for each sample
% Row 4 = List of indentation numbers in a group including edited ones




function final_modulus = auto_accept_curve(path_name,filenum,indenter,C_S)
global min_depth max_depth stp_depth I_R

if filenum > 99
    namedeguer = num2str(filenum);
elseif (filenum >= 10) && (filenum < 100)
    namedeguer = ['0' num2str(filenum)];
elseif filenum < 10
    namedeguer = ['00' num2str(filenum)];
else
    error('File number must be between 0-999');
    return
end
file_list = dir([path_name '*' namedeguer '*.fv001lq']);
% Kill program if the directory listing returns an empty handle
if isempty(file_list) == 1
    error('Requested file(s) cannot be found!');
    return
end

% Indenter types are
%       1-standard blunt cone
%       2-custom sphere
%       3-custom hemisphere
if indenter == 1
    R_tip = 50e-9; % Radius of curvature
    t_tip = 40*pi/180; % Cone semi-angle
    b_tip = R_tip*cos(t_tip); % Cylindrical radius
elseif indenter == 2
    % Indenter radius
elseif indenter == 3
    % Hemisphere radius
end

% Import the analysis result file form the SGI:
rawdata = dlmread([path_name file_list(1).name], '\t', 3, 0);
% Create an initial matrix of extension-deflection vectors:
def_ex = rawdata(:,1:2);
% Import the contact point computed by the SGI:
contact_point = rawdata(7,21);
% Parse the initial vectors into pre- and post-contact:
index_contact = find(def_ex(:,1) >= contact_point,1);
pre_contact = def_ex(1:index_contact,1:2);
post_contact = def_ex(index_contact:end,1:2);

% FORCE-DEPTH CALCULATION
% Calculate force-depth relationship based on the input values for
% the cantilever stiffness:
force = (post_contact(:,2) - post_contact(1,2)) .* 1e-9 * C_S;
depth = ((post_contact(:,1) - post_contact(1,1)) .* 1e-6) ...
    - ((post_contact(:,2) - post_contact(1,2)) .* 1e-9);

% POINTWISE MODULUS CALCULATION
% Calculate depth-dependent pointwise modulus based on the input
% values for the indenter type and geometry:
if indenter == 1
    for m = 1:length(depth)
        h_tip = depth(m);
        a_tip = sqrt(R_tip * h_tip);
        % The solution for the modulus uses a spherical geometry
        % until the contact radius surpasses the cylindrical radius
        % (b_tip) and a blunt-cone geometry once this point is
        % reached:
        if a_tip <= b_tip
            modulus(m,1) = (0.75*force(m))/...
                (sqrt(R_tip)*(h_tip^1.5));
        elseif a_tip > b_tip
            % Summon the GET_CONTACT_RADIUS function that
            % numerically solves for the contact radius using
            % equation (4) of Briscoe et al. 1999:
            a_tip = get_contact_radius(b_tip,h_tip,R_tip,t_tip);
            messy_part = (a_tip*h_tip - (a_tip^2/(2*tan(t_tip)))*...
                (0.5*pi - asin(b_tip/a_tip)) - (a_tip^3/(3*R_tip))...
                + sqrt(a_tip^2 - b_tip^2)*((b_tip/(2*tan(t_tip)))...
                + (a_tip^2 - b_tip^2)/(3*R_tip)));
            modulus(m,1) = (0.5*force(m))/messy_part;
        end
        a_vector(m,1) = a_tip;
    end
    fit_E = NaN;
elseif indenter == 2
    modulus = (3/4) .* force .* (I_R^-0.5) .* (depth.^-1.5);
    for dp = 1:1:length(depth)
        if depth(dp) <= 0.9*max(depth)
            tobefit_depth(dp,1) = depth(dp);
            tobefit_force(dp,1) = force(dp);
        end
    end
    fit_E = ((I_R^0.5) .* (tobefit_depth.^1.5)) \ ...
        ((3/4) .* tobefit_force);
    fit_force = (4/3) * fit_E .* (I_R^0.5) .* (tobefit_depth.^1.5);
elseif indenter == 3
    modulus = force .* I_R^-1;
    fit_E = I_R \ tobefit_force;
    fit_force = fit_E .* I_R;
end


depth = depth*1e9; m = 1; % Current row number
for k = min_depth:stp_depth:max_depth
    e_count = ((k-min_depth)/stp_depth) + 1;
    if k < max(depth)
        % Increment the rows until the next row reads a value greater
        % than the requested depth (and the current number reads less
        % than or equal to the requested depth)
        while depth(m+1) < k
            m = m + 1;
        end
        % If the condition is met interpolate for the pointwise modulus
        if (depth(m) < k) && depth(m+1) >= k
            % Below is the interpolation function that finds the exact
            % pointwise modulus at this given depth:
            %        y_int = y1 + (x_int - x1)*((y2-y1)/(x2-x1))
            final_modulus(e_count,1) = modulus(m) + (k - depth(m))*...
                ((modulus(m+1)-modulus(m))/(depth(m+1)-depth(m)));
        end     
    else
        final_modulus(e_count,1) = NaN;
    end     
end
final_modulus(e_count+1,1) = fit_E;





function final_modulus = user_accept_curve(path_name,filenum,indenter,C_S)
global min_depth max_depth stp_depth I_R

if filenum > 99
    namedeguer = num2str(filenum);
elseif (filenum >= 10) && (filenum < 100)
    namedeguer = ['0' num2str(filenum)];
elseif filenum < 10
    namedeguer = ['00' num2str(filenum)];
else
    error('File number must be between 0-999');
    return
end
file_list = dir([path_name '*' namedeguer '*.fv001lq']);
% Kill program if the directory listing returns an empty handle
if isempty(file_list) == 1
    error('Requested file(s) cannot be found!');
    return
end

% Indenter types are
%       1-standard blunt cone
%       2-custom sphere
%       3-custom hemisphere
if indenter == 1
    R_tip = 50e-9; % Radius of curvature
    t_tip = 40*pi/180; % Cone semi-angle
    b_tip = R_tip*cos(t_tip); % Cylindrical radius
elseif indenter == 2
    % Indenter radius
elseif indenter == 3
    % Hemisphere radius
end


% Figure position coordinates to be used later:
% ----- (from B.S. Elkin)
% create position vectors to position the two figures in the top
% third of the screen 
bdwidth = 5;
topbdwidth = 30;
set(0,'Units','pixels') 
scnsize = get(0,'ScreenSize');
pos1  = [bdwidth, 1/2*scnsize(4) + bdwidth, ...
    scnsize(3)/2 - 2*bdwidth, scnsize(4)/2 - ...
    (topbdwidth + bdwidth)];
pos2 = [pos1(1) + scnsize(3)/2, pos1(2), pos1(3), pos1(4)];
% ----- (end quote)


check_opt = 2;
while (check_opt == 2) || (check_opt == -1)
    % Import the analysis result file form the SGI:
    rawdata = dlmread([path_name file_list(1).name], '\t', 3, 0);
    % Create an initial matrix of extension-deflection vectors:
    def_ex = rawdata(:,1:2);
    % Import the contact point computed by the SGI:
    costa_contact_point = rawdata(7,21);
    if check_opt == 2
        contact_point = costa_contact_point;
    else
        contact_point = hand_selected_point;
    end
    % Parse the initial vectors into pre- and post-contact:
    index_contact = find(def_ex(:,1) >= contact_point,1);
    pre_contact = def_ex(1:index_contact,1:2);
    post_contact = def_ex(index_contact:end,1:2);

    % FORCE-DEPTH CALCULATION
    % Calculate force-depth relationship based on the input values for
    % the cantilever stiffness:
    force = (post_contact(:,2) - post_contact(1,2)) .* 1e-9 * C_S;
    depth = ((post_contact(:,1) - post_contact(1,1)) .* 1e-6) ...
        - ((post_contact(:,2) - post_contact(1,2)) .* 1e-9);

    % POINTWISE MODULUS CALCULATION
    % Calculate depth-dependent pointwise modulus based on the input
    % values for the indenter type and geometry:
    if indenter == 1
        for m = 1:length(depth)
            h_tip = depth(m);
            a_tip = sqrt(R_tip * h_tip);
            % The solution for the modulus uses a spherical geometry
            % until the contact radius surpasses the cylindrical radius
            % (b_tip) and a blunt-cone geometry once this point is
            % reached:
            if a_tip <= b_tip
                modulus(m,1) = (0.75*force(m))/...
                    (sqrt(R_tip)*(h_tip^1.5));
            elseif a_tip > b_tip
                % Summon the GET_CONTACT_RADIUS function that
                % numerically solves for the contact radius using
                % equation (4) of Briscoe et al. 1999:
                a_tip = get_contact_radius(b_tip,h_tip,R_tip,t_tip);
                messy_part = (a_tip*h_tip - (a_tip^2/(2*tan(t_tip)))*...
                    (0.5*pi - asin(b_tip/a_tip)) - (a_tip^3/(3*R_tip))...
                    + sqrt(a_tip^2 - b_tip^2)*((b_tip/(2*tan(t_tip)))...
                    + (a_tip^2 - b_tip^2)/(3*R_tip)));
                modulus(m,1) = (0.5*force(m))/messy_part;
            end
            a_vector(m,1) = a_tip;
        end                
        fit_E = NaN;
    elseif indenter == 2
        modulus = (3/4) .* force .* (I_R^-0.5) .* (depth.^-1.5);
        for dp = 1:1:length(depth)
            if depth(dp) <= 0.9*max(depth)
                tobefit_depth(dp,1) = depth(dp);
                tobefit_force(dp,1) = force(dp);
            end
        end
        fit_E = ((I_R^0.5) .* (tobefit_depth.^1.5)) \ ...
            ((3/4) .* tobefit_force);
        fit_force = (4/3) * fit_E .* (I_R^0.5) .* (tobefit_depth.^1.5);
    elseif indenter == 3
        modulus = force .* I_R^-1;
        fit_E = I_R \ tobefit_force;
        fit_force = fit_E .* I_R;
    end

    % Plot the processed data including the raw force curve, denoting the
    % contact point by a color change, the fitted force curve, and modulus
    % + force vs. depth for pointwise analysis.
    %
    % Important Note: Currently, the conic indenter is not fitted to give a
    % single modulus value, therefore a fit result cannot be plotted for
    % that scenario. The below if-statement is written for that purpose:
    if indenter == 1
        f_plot = figure('Position',pos1);
        plot(pre_contact(:,1),pre_contact(:,2),'-b',...
            post_contact(:,1),post_contact(:,2),'-r');
        title(['Raw Force Curve: ' file_list(1).name '']);
        m_plot = figure('Position',pos2);
        subplot(2,1,1),plot((depth.*1e9),(force.*1e9)); 
        ylabel('Force (nN)'),xlabel('Depth (nm)');
        subplot(2,1,2),plot((depth.*1e9),(modulus.*1e-3));
        ylabel('Modulus (kPa)'),xlabel('Depth (nm)');
    else
         f_plot = figure('Position',pos1);
        plot(pre_contact(:,1),pre_contact(:,2),'-b',...
            post_contact(:,1),post_contact(:,2),'-r');
        title(['Raw Force Curve: ' file_list(1).name '']);
        m_plot = figure('Position',pos2);
        subplot(2,1,1),plot((depth.*1e9),(force.*1e9),'-b',...
            (tobefit_depth.*1e9),(fit_force.*1e9),'--g'); 
        ylabel('Force (nN)'),xlabel('Depth (nm)');
        subplot(2,1,2),plot((depth.*1e9),(modulus.*1e-3));
        ylabel('Modulus (kPa)'),xlabel('Depth (nm)');
    end
    
    % Ask the user to confirm the contact point:
    check_button = questdlg('Do you want to accept force data?',...
        'Edit Contact Point','Yes','No','Edit','Yes');
    if strcmp(check_button,'Yes')
        check_opt = 0;
    elseif strcmp(check_button,'No')
        check_opt = 1;
        modulus(:,1) = NaN;
    elseif strcmp(check_button,'Edit')
        check_opt = -1;
    end
    clear check_button

    if check_opt == -1
       figure(f_plot)
       [hand_selected_point,temp] = ginput(1);
       clear temp depth modulus force
       close(f_plot)
       close(m_plot)
   end
end

close(f_plot)
close(m_plot)

depth = depth*1e9; m = 1; % Current row number
for k = min_depth:stp_depth:max_depth
    e_count = ((k-min_depth)/stp_depth) + 1;
    if k < max(depth)
        % Increment the rows until the next row reads a value greater
        % than the requested depth (and the current number reads less
        % than or equal to the requested depth)
        while depth(m+1) < k
            m = m + 1;
        end
        % If the condition is met interpolate for the pointwise modulus
        if (depth(m) < k) && depth(m+1) >= k
            % Below is the interpolation function that finds the exact
            % pointwise modulus at this given depth:
            %        y_int = y1 + (x_int - x1)*((y2-y1)/(x2-x1))
            final_modulus(e_count,1) = modulus(m) + (k - depth(m))*...
                ((modulus(m+1)-modulus(m))/(depth(m+1)-depth(m)));
        end     
    else
        final_modulus(e_count,1) = NaN;
    end     
end
final_modulus(e_count+1,1) = fit_E;
