function [curves,total_error] = afm_post_pro(path_name,file_name,...
    geometry_opt,pointwise_opt)
% function [curves,total_error] = afm_post_pro(geometry_opt,pointwise_opt)
%   Produces a CURVES matrix that is composed of 7 key components:
%   1. X-Coordinate of the indentation
%   2. Y-Coordinate of the indentation
%   3. Extension curve (distance-force in 'm' and 'N')
%   4. Retraction curve (distance-force in 'm' and'N')
%   5. Force curve (depth-force in 'nm' and 'nN')
%   6. Pointwise modulus curve (depth-modulus in 'nm' and 'kPa')
%   7. Fitted apparent modulus (kPa)
%
% AZELOGLU E.U. (c) 2016




% =========================================================================
% clear;
if nargin < 2
    [file_name,path_name] =  uigetfile({'*.*','All Files (*.*)',...
        }, 'Select an indentation file to view');
    if file_name == 0
        clear file_name
        clear path_name
        return
    end
else
    path_name = [path_name '/'];
end
fid=fopen([path_name,file_name],'r');
tLines = fgets(fid);
% numCols = numel(strfind(tLines,',')) + 1;
file_header = strsplit(tLines,',');
no_of_curves = (length(file_header) - 2)/6;
no_of_curves = 16; % Temporary
total_error = cell(no_of_curves,1);

for n = 1:no_of_curves
    column_no = 6*(n-1) + 1;
    line_no_temp = strsplit(file_header{1,column_no},'Point');
    curves{n,1} = str2num(strrep(line_no_temp{1,1},'Line','')) + 1;
    curves{n,2} = str2num(strrep(line_no_temp{1,2},'Ind_Ext','')) + 1;
    fmt_temp = [ repmat('%*f',1,column_no-1), '%f%f%f%f%[^\n]'];
    fid=fopen([path_name,file_name],'r'); 
    data_temp = textscan(fid,fmt_temp,'Delimiter',',','HeaderLines',1);
    ext_temp = [data_temp{1,1} data_temp{1,2}];
    ret_temp = [data_temp{1,3} data_temp{1,4}];
    curves{n,3} = ext_temp;
    curves{n,4} = ret_temp;
    clear *_temp;
    fclose(fid);
end


% -------------------------------------------------------------------------




% =========================================================================
% PROGRAM OPTIONS
% If the deflection sensitivity recorded in the file is wrong or needs to
% be corrected on the spot then the user can interfere by plugging in the
% correct deflection sensitivity value:
defl_sens_opt = 0;
% If defl_sens_opt is ZERO then the code uses the encoded deflection
% sensitivity for the file, if the option is set to ONE then the code
% prompts the user to input the correct deflection sensitivity.

% Load-Unload option chooses which portion of the force curve to analyze.
% If the option is set to ONE, then advance is loaded into the force curve;
% if the option is set to TWO, the retract portion is loaded:
extension_opt = 1; % Make this a variable

% Number of points to be skipped for contact point determination:
NSKIP = 100; % Make this a variable

% The depth of data to be included in the contact point fitting:
% The y_crit is in N. It's empirically selected as 1nN, which is
% approximately equivalent to a 10nm deflection on a 0.1nN/nm cantilever
% that is commonly used for cell indentation.
y_crit = 1e-9;  % Make this a variable

% Orders of the pre-contact and post-contact fits:
N_post = 2; % Make this a variable
N_pre = 1; % Make this a variable

% To plot the individual graphs and create a movie of the process set
% movie_opt to 1, otherwise set it to 0:
movie_opt = 0; % Make this a variable
% Name of the movie file:
if movie_opt == 1
    aviobj = avifile('error_pro_new.avi','quality',90); % Make a variable
end

% To equate the two functions at the contact point, set equate_opt to 1, if
% the analytical solution for linear-quadratic fit is to be used then set
% it to 2, otherwise set it to 0:
equate_opt = 1; % Make this a variable

% Real-time computation display option:
display_opt = 0; % Make this a variable
% If display_opt is zero, then the code will just show the current step; if
% the option is set to one, then it will show the current step of
% calculation (n_current / n_total); finally, if display_opt is set to two,
% the code will display the depth it is currently analyzing.

% An additional mode where the code will try the neighboring points of the
% optimized contact point with equal lengths of pre- and post-contact
% curves to better see the effects of contact point assumption:
shift_test_mode = 0; % Make this a variable

% Cantilever stiffness in N/m:
cantilever_stiffness = 0.1; % Make this a variable

% Tip geometry selection option:
geometry_opt = 1; % Make this a variable
% Geometry option determines the kind of contact geometry that is used to
% calculate the modulus. If geometry_opt is ONE, then the standard Si3N4
% tip is chosen; if the option is set to TWO, then custom spherical probe
% geometry is assumed; if the option is set to THREE, then a flat
% cylindrical (hemispherical) tip geometry is assumed; finally, if the
% option is set to FOUR the generalized contact geometry case is utilized
% with completely user defined parameters.

% Depthwise modulus calculation option:
pointwise_opt = 1; % Make this a variable
% If pointwise_opt is zero, then the code ignores the pointwise calculation
% kernel; if the option is set to one, then it calculates the modulus only
% at the marked locations (multiples of fifty nanometers of depth);
% finally, if the option is set to two, the code calculates modulus at each
% location.
% -------------------------------------------------------------------------



% =========================================================================
clc; disp('Pre-processing data...');

for c_n = 1:no_of_curves
    disp(['Processing ' num2str(c_n) ' out of ' ... 
        num2str(no_of_curves) '.']);
    extension = [curves{c_n,3}];
    retract = [curves{c_n,4}];
    % ---------------------------------------------------------------------



    % =====================================================================
    % CONTACT POINT DETERMINATION
    % Use Nth order polynomials, determined by the user, to fit the pre- and
    % post-contact regions respectively, in order to detect the contact point.
    % The algorithm shall skip the first "NSKIP" number of data points (again
    % determined by the user) and then create a new data set from this point
    % until it reaches a user defined "y_crit" amount of deflection. The new
    % matrix will be recorded from the selected portion of the "extension"
    % matrix until "y_not + y_crit" is reached or the end-of-file.
    y_not = mean(extension(NSKIP:NSKIP+10,2));


    % ---------------------------------------------------------------------
    % Start the post-contact fitting process from 50% of the maximum
    % deflection threshold line and progress towards the pre-contact until the
    % sum of residuals of the two polynomial fits are minimized. Designate 5%
    % of the maximum deflection, or 50% increase in error - whichever comes
    % first - as stop markers for the fitting process. Record the sum of the
    % residuals to a comparison vector and find the location of the minimum on
    % that vector.
    end_mark = 0; beg_mark = 0;
    if max(extension(:,2)) <= y_not + y_crit
        warning(['Not enough post-contact. Total indent force < ' ...
            num2str(y_crit) 'N']);
    elseif isnan(extension(NSKIP,2)) == 1
        warning('Not enough pre-contact.');
    else
        if extension_opt == 1
            samps_line = length(extension);
        elseif extension_opt == 2
            samps_line = length(retract);
        end

        for n_fit = [NSKIP:1:samps_line]
            if extension_opt == 1
                fit_to_this(n_fit-NSKIP+1,1) = extension(n_fit+1,1);
                fit_to_this(n_fit-NSKIP+1,2) = extension(n_fit+1,2);
            elseif extension_opt == 2
                fit_to_this(n_fit-NSKIP+1,1) = retract(n_fit+1,1);
                fit_to_this(n_fit-NSKIP+1,2) = retract(n_fit+1,2);
            end

            if (fit_to_this(n_fit-NSKIP+1,2) >= y_not + (0.05*y_crit)) &&...
                    isequal(end_mark,0) == 1
                end_mark = n_fit - NSKIP + 1;
            elseif (fit_to_this(n_fit-NSKIP+1,2) >= y_not + (0.3*y_crit)) &&...
                    isequal(beg_mark,0) == 1
                beg_mark = n_fit - NSKIP + 1;
            elseif fit_to_this(n_fit-NSKIP+1,2) >= y_not + y_crit
                fit_size = n_fit - NSKIP + 1;
                break
            end
        end
        clear n_fit

        % ---------------------------------------------------------------------
        disp('Calculating optimum contact point...');
        % START the FITTING:
        for n = [beg_mark:-1:end_mark]
            m = (beg_mark - n) + 1;
            if display_opt == 1
                disp(['Fitting... ' num2str(m) '/' ...
                    num2str(beg_mark-end_mark+1)]);
            elseif display_opt == 2
                disp(['Fitting... ' num2str(fit_to_this(n,1)) ' um']);
            end

            % Divide the selected portion vector, "fit_to_this", into pre- and
            % post-contact portions:
            pre_fit = fit_to_this(1:n,:);
            post_fit = fit_to_this(n:fit_size,:);

            % Create a depth vector to keep track of location:
            z_fit_time(m) = fit_to_this(n,1);

            % Create a coordinate vector to determine the final coordinate for
            % the x* value with minimum error:
            x_star_coord(m) = n;

                % -------------------------------------------------------------
                % DO NOT EQUATE FUNCTIONS AT X*:
                % The first loop does not equate the two functions at the 
                % contact point. It minimizes the errors independently.
            if equate_opt == 0
                % Fitting for the pre-contact portion:
                pre_poly = polyfit(pre_fit(:,1),pre_fit(:,2),N_pre);
                pre_f = polyval(pre_poly,pre_fit(:,1))';
                pre_err(m) = sum(sqrt((pre_f(:) - pre_fit(:,2)).^2));
                % Fittting for the post-contact portion:
                post_poly = polyfit(post_fit(:,1),post_fit(:,2),N_post);
                post_f = polyval(post_poly,post_fit(:,1))';
                post_err(m) = sum(sqrt((post_f(:) - post_fit(:,2)).^2));

                total_err = pre_err + post_err;

                % The loop below, if opted, will generate a movie in the 'AFM 
                % Post Processor' directory under the designated name.
                if movie_opt == 1
                    figure;
                    plot(pre_fit(:,1),pre_fit(:,2),'.w',pre_fit(:,1),pre_f,'-w',...
                        post_fit(:,1),post_fit(:,2),'.y',post_fit(:,1),post_f,'-y');
                    str1(1) = {['Contact Point: ' num2str(fit_to_this(n,1)) ...
                        ' | Pre-Error: '...
                            num2str(pre_err(n-end_mark+1)) ' | Post-Error: '...
                            num2str(post_err(n-end_mark+1)) '']};
                    str1(2) = {['Total Error: ' num2str(pre_err(n-end_mark+1)+...
                                post_err(n-end_mark+1)) '']};
                    uicontrol('Style','text','BackgroundColor','k',...
                        'ForegroundColor','w',...
                        'Position',[132 60 350 32],'String',str1);
                    title('Live Fitting Progress Viewer');
                    F = getframe(gca);
                    aviobj = addframe(aviobj,F);
                end

            % -----------------------------------------------------------------
            % EQUATE FUNCTIONS AT X*
            % The second set of loops are used if the pre- and post-contact curves
            % are forced to be equal at the chosen contact point. These loops go
            % through a very different optimization process.
            elseif equate_opt == 1
                % The first way of fitting the two functions uses the constrained
                % minimization tool that minimizes the RMS error between the
                % fitted values and the data point while keeping the contact points
                % of the pre-contact and post-contact portions equal. This portion
                % of the code summons the "constrained_dual_fit" function that uses
                % the built-in MATLAB function "fmincon" to find the coefficients
                % while using the unconstraint numbers as startup guess values.

                [fit_coefs,total_rms,convergence_flag] = ...
                    constrained_dual_fit(pre_fit,post_fit,N_pre,N_post);

                if convergence_flag < 0
                    warning('The fit did not converge to a solution')
                end

                % Load the results for the pre-contact curve into a cell array:
                pre_f{1,m} = pre_fit(:,1);
                pre_f{2,m} = pre_fit(:,2);
                % Evaluate the pre-contact function based on the power of fits
                % entered by the user:
                if N_pre == 2
                    pre_f{3,m} = polyval(...
                        [fit_coefs(1) fit_coefs(2) fit_coefs(3)],pre_fit(:,1));
                elseif N_pre == 1
                    pre_f{3,m} = polyval([fit_coefs(1) fit_coefs(2)],...
                        pre_fit(:,1));
                end
                pre_err(m) = sum((pre_fit(:,2)-pre_f{3,m}).^2);

                % Load the results for the post-contact curve into a cell array:
                post_f{1,m} = post_fit(:,1);
                post_f{2,m} = post_fit(:,2);
                if N_post == 3
                    post_f{3,m} = polyval([fit_coefs(N_pre+2) fit_coefs(N_pre+3)...
                            fit_coefs(N_pre+4) fit_coefs(N_pre+5)],...
                        (post_fit(:,1)-post_fit(1,1)));
                elseif N_post == 2
                    post_f{3,m} = polyval([fit_coefs(N_pre+2) fit_coefs(N_pre+3)...
                            fit_coefs(N_pre+4)],(post_fit(:,1)-post_fit(1,1)));
                elseif N_post == 1
                    post_f{3,m} = polyval(...
                        [fit_coefs(N_pre+2) fit_coefs(N_pre+3)], ...
                        (post_fit(:,1)-post_fit(1,1)));
                end
                post_err(m) = sum((post_fit(:,2)-post_f{3,m}).^2);

                total_err(m) = total_rms;

                % The loop below, if opted, will generate a movie in the 'AFM
                % Post Processor' directory under the designated name.
                if movie_opt == 1
                    figure;
                    plot(pre_fit(:,1),pre_fit(:,2),'.w',...
                        pre_fit(:,1),pre_f{3,m},'-w',...
                        post_fit(:,1),post_fit(:,2),'.y',...
                        post_fit(:,1),post_f{3,m},'-y');
                    str1(1) = {['Contact Point: ' num2str(fit_to_this(n,1)) '']};
                    str1(2) = {['Total Error: ' num2str(total_err(m)) '']};
                    uicontrol('Style','text','BackgroundColor','k',...
                        'ForegroundColor','w',...
                        'Position',[132 60 350 32],'String',str1);
                    title('Live Fitting Progress Viewer');
                    F = getframe(gca);
                    aviobj = addframe(aviobj,F);
                end

                if total_err(m) <= min(total_err)
                    min_err = [total_err(m) m];
                    if min_err(1) > 1e6
                        disp('Fit does not converge to a solution');
                        disp('Switching to the analytical lin-quad solution...');
                        equate_opt = 2;
                    end
                elseif (total_err(m) >= (1.33 * min_err(1))) && (m > 10)
                    disp(['Minimum error achieved: ' num2str(min_err(1)) '']);
                    break
                end
            end


            % -----------------------------------------------------------------
            % EQUATE FUNCTIONS AT X* (METHOD #2)     
            if equate_opt == 2
                % The code below will simultaneously fit functions of first and
                % second orders to pre- and post-contact sections, respectively.

                % The fits will be forced to be equal at an optimal height (y*) at
                % the prescribed contact point (x*). The protocol is carried out by
                % calling the "fit_pre_post" function that uses a least-squares
                % approach to analytically solve the over-determined constraint
                % optimization problem. This is different than the previous
                % approach since it uses an analytical solution instead of an
                % optimization solution. This approach has a stronger chance of
                % convergence; therefore it's used as a back-up module.

                % Summon FIT_PRE_POST function to fit the data:
                [fit_coefs,fit_rms] = fit_pre_post(pre_fit,post_fit,...
                    N_pre,N_post);
                % Load the results for the pre-contact curve into a cell array:
                pre_f{1,m} = pre_fit(:,1);
                pre_f{2,m} = pre_fit(:,2);
                pre_f{3,m} = polyval([fit_coefs(2) fit_coefs(1)],pre_fit(:,1));
                pre_err(m) = sum((pre_fit(:,2)-pre_f{3,m}).^2);
                % Load the results for the post-contact curve into a cell array:
                post_f{1,m} = post_fit(:,1);
                post_f{2,m} = post_fit(:,2);
                post_f{3,m} = polyval([fit_coefs(5) fit_coefs(4) fit_coefs(3)],...
                    (post_fit(:,1)-post_fit(1,1)));
                post_err(m) = sum((post_fit(:,2)-post_f{3,m}).^2);

                total_err = pre_err + post_err;

                % The loop below, if opted, will generate a movie in the 'AFM
                % Post Processor' directory under the designated name.
                if movie_opt == 1
                    figure;
                    plot(pre_fit(:,1),pre_fit(:,2),'.w',...
                        pre_fit(:,1),pre_f{3,m},'-w',...
                        post_fit(:,1),post_fit(:,2),'.y',...
                        post_fit(:,1),post_f{3,m},'-y');
                    str1(1) = {['Contact Point: ' num2str(fit_to_this(n,1)) '']};
                    str1(2) = {['Total Error: ' num2str(total_err(m)) '']};
                    uicontrol('Style','text','BackgroundColor','k',...
                        'ForegroundColor','w',...
                        'Position',[132 60 350 32],'String',str1);
                    title('Live Fitting Progress Viewer');
                    F = getframe(gca);
                    aviobj = addframe(aviobj,F);
                end

                if total_err(m) <= min(total_err)
                    min_err = [total_err(m) m];
                    if min_err(1) > 1e6
                        error('Fit does not converge to a solution');
                    end
                elseif (total_err(m) >= (1.33 * min_err(1))) && (m > 10)
                    disp(['Minimum error achieved: ' num2str(min_err(1)) '']);
                    break
                end
            end
        end

        if movie_opt == 1
            aviobj = close(aviobj);
        end

        [C I] = min(total_err);
        total_err = [z_fit_time' total_err'];
        if extension_opt == 1
            pre_final = extension(1:(NSKIP+x_star_coord(I)),:);
            post_final = extension((NSKIP+x_star_coord(I)+1):end,:);
        elseif extension_opt == 2
            pre_final = retract(1:(NSKIP+x_star_coord(I)),:);
            post_final = retract((NSKIP+x_star_coord(I)+1):end,:);
        end
        clear z_fit_time fit_*
        % ---------------------------------------------------------------------



        % =====================================================================
        % Create a FORCE-DEPTH CURVE (FC):

        % Find the deflection and force from the post-contact curve and load them
        % into the FC matrix. Since the units we want to use in the FC are in
        % nanometers, we convert the micron units of the post-contact curves here.
        FC = zeros(size(post_final,1),2);
        for n=1:size(post_final,1)
            FC(n,1) = (post_final(n,1) - post_final(1,1))*1e9;
            if post_final(1,2) < 0
                FC(n,2) = 1e9*(post_final(n,2)-post_final(1,2));
            else
                FC(n,2) = 1e9*post_final(n,2);
            end
        end
        % The final form of the force curve has the depth in its first column and
        % force in its second column. The third and fourth columns are
        % inconsequential.

        if geometry_opt == 1
            tip_geo = 'Standard Si3N4';
        elseif geometry_opt == 2
            tip_geo = 'Spherical';
        elseif geometry_opt == 3
            tip_geo = 'Cylindrical';
        else
            tip_geo = 'Custom';
        end
        % ---------------------------------------------------------------------



        % =====================================================================
        % APPARENT MODULUS CALCULATION
        disp('Calculating apparent modulus and pointwise moduli...');
        % Summon the "calculate_moduli" function to compute elastic modulus:
        if pointwise_opt == 0
            E_apparent = ...
                calculate_moduli(FC,geometry_opt,pointwise_opt,display_opt);
        else
            [E_apparent,E_pointwise] = ...
                calculate_moduli(FC,geometry_opt,pointwise_opt,display_opt);
        end
        % ---------------------------------------------------------------------

        colordef black;


        % =====================================================================
        % TEST PROXIMITY CONTACT POINTS
        if shift_test_mode == 1
            figure;
            pre_post_length = floor(samps_line/10);
            width_to_shift = 20; 
            k = 1;
            for n_temp = 0:(2*width_to_shift)
                move_x_star = n_temp - width_to_shift;
                if move_x_star >= -width_to_shift && move_x_star < -width_to_shift/2
                    colour = '-y'; quarter = 1; c_colour = '+y';
                elseif move_x_star >= -width_to_shift && move_x_star < 0
                    colour = '-r'; quarter = 2; c_colour = '+r';
                elseif move_x_star >= 0 && move_x_star < width_to_shift/2
                    colour = '-c'; quarter = 3; c_colour = '+c';
                elseif move_x_star > width_to_shift/2
                    colour = '-g'; quarter = 4; c_colour = '+g';
                end

                post_start = NSKIP+x_star_coord(I)+move_x_star+1;
                test_post = extension(post_start:(post_start+pre_post_length),1:2);
                pre_start = NSKIP+x_star_coord(I)+move_x_star;
                test_pre = extension((pre_start-pre_post_length):pre_start,1:2);
                [temp_coefs,temp_rms,temp_flag] = ...
                    constrained_dual_fit(test_pre,test_post,N_pre,N_post);

                if N_pre == 2
                    test_pre_fit = polyval([temp_coefs(1) ...
                        temp_coefs(2) temp_coefs(3)], test_pre);
                elseif N_pre == 1
                    test_pre_fit = polyval([temp_coefs(1) temp_coefs(2)],...
                        test_pre(:,1));
                end
                if N_post == 3
                    test_post_fit = polyval(...
                        [temp_coefs(N_pre+2) temp_coefs(N_pre+3) ...
                        temp_coefs(N_pre+4) temp_coefs(N_pre+5)],...
                        (test_post(:,1)-test_post(1,1)));
                elseif N_post == 2
                    test_post_fit = polyval(...
                        [temp_coefs(N_pre+2) temp_coefs(N_pre+3) ...
                        temp_coefs(N_pre+4)],(test_post(:,1)-test_post(1,1)));
                elseif N_post == 1
                    test_post_fit = polyval(...
                        [temp_coefs(N_pre+2) temp_coefs(N_pre+3)], ...
                        (test_post(:,1)-test_post(1,1)));
                end

                if k == 1
                    subplot(1,4,1),plot(test_pre(:,1),...
                        test_pre(:,2),'.w',test_post(:,1),test_post(:,2),...
                        '.m',test_post(1,1),test_post(1,2),c_colour); hold on
                else
                    subplot(1,4,1),plot(test_post(1,1),test_post(1,2),c_colour); 
                    hold on
                end
        %          subplot(1,4,1),plot(test_pre(:,1),...
        %              test_pre(:,2),'.w',test_pre(:,1),...
        %              test_pre_fit,colour,test_post(:,1),test_post(:,2),'.w',...
        %              test_post(:,1),test_post_fit,colour);
                set(gca,'XGrid','on','YGrid','on');

                end_post_test = extension(post_start:end,:);
                end_samp_no = size(end_post_test,1);
                for n=1:end_samp_no
                    FC_temp(n,3) = 1000 * end_post_test(end+1-n,1);
                    if extension_opt == 1
                        FC_temp(n,4) = end_post_test(n,2) - end_post_test(1,2);
                    elseif extension_opt == 2
                        FC_temp(n,4) = end_post_test(n,3) - end_post_test(1,3);
                    end
                    FC_temp(n,1) = FC_temp(n,3) - FC_temp(n,4);
                    FC_temp(n,2) = FC_temp(n,4) * cantilever_stiffness;
                end

                if pointwise_opt == 0
                    E_app_temp = ...
                        calculate_moduli(FC_temp,geometry_opt,...
                        pointwise_opt,display_opt);
                else
                    [E_app_temp,E_p_temp] = ...
                        calculate_moduli(FC_temp,geometry_opt,...
                        pointwise_opt,display_opt);
                end
                E_variation = 100 * (std(E_p_temp(2:end,2)) /...
                    mean(E_p_temp(2:end,2)));
                E_slope = polyfit(E_p_temp(4:end,1),E_p_temp(4:end,2),1);
                E_slope = E_slope * 1000;
                subplot(1,4,2),plot(FC_temp(:,1),FC_temp(:,2),colour); 
                set(gca,'XGrid','on','YGrid','on'); hold on;
                subplot(1,4,3),plot(E_p_temp(:,1),E_p_temp(:,2),colour);
                set(gca,'XGrid','on','YGrid','on'); hold on;

                if rem(k,4) == 0
                    subplot(1,4,4); axis off;
                    information = {['x*: ' num2str(move_x_star)] ...
                            %['Modulus: ' num2str(E_app_temp) ' KPa'] ...
                            %['Deviation from Mean E: ' num2str(E_variation) '%']; ...
                            ['Slope of E_a_p_p(D): ' num2str(E_slope(1))]};
                            %['Minimum Error: ' num2str(temp_rms)]};
                    text(0.05,0.025*k,information);
                end

                proximity_results{k,1} = move_x_star;
                proximity_results{k,2} = temp_coefs;
                proximity_results{k,3} = temp_rms;
                proximity_results{k,4} = E_slope(1);
                proximity_results{k,5} = E_app_temp;
                proximity_results{k,6} = E_p_temp;
                proximity_results{k,7} = FC_temp;
                clear E_p_temp FC_temp
                k = k + 1;
            end
            set(gcf,'Position',[30 60 930 500]);
        end

        % =====================================================================
        % PLOTTING THE RESULTS
        if shift_test_mode == 0      
            figure;
            
            subplot(2,3,1),plot(extension(:,1),extension(:,2),... % Extend
                '-w', retract(:,1),retract(:,2),'-y');            % Retract
            title('Raw AFM Deflection Curve');
            set(gca,'XGrid','on','YGrid','on');
            desired_ax = axis;

            if extension_opt == 1
                subplot(2,3,2),plot(pre_final(:,1),pre_final(:,2),'.c',...
                    post_final(:,1),post_final(:,2),'.m');
            elseif extension_opt == 2
                subplot(2,3,2),plot(pre_final(:,1),pre_final(:,3),'.c',...
                    post_final(:,1),post_final(:,3),'.m');
            end
            title('Optimum Contact Point');
            legend('Pre-Contact Curve','Post-Contact Curve');
            set(gca,'XGrid','on','YGrid','on');
            axis(desired_ax);

            subplot(2,3,3),plot(total_err(:,1),total_err(:,2),'.w');
            title('Total RMS Error for Prescribed X*');
            set(gca,'XGrid','on','YGrid','on','XLim',[desired_ax(1:2)]);

            subplot(2,3,4),plot(FC(:,1),FC(:,2),'-w');
            title('Force vs. Indentation Depth');
            set(gca,'XGrid','on','YGrid','on');

            subplot(2,3,5),plot(E_pointwise(:,1),E_pointwise(:,2),'-w');
            title('Pointwise Apparent Modulus vs. Indentation Depth');
            set(gca,'XGrid','on','YGrid','on');

            subplot(2,3,6); axis off;
            information = {['File name: ' file_name ]; ...
                    ['Tip geometry: ' tip_geo ]; ...
                    ['Fitted modulus: ' num2str(E_apparent) ' KPa']; ...
                    ['Minimum Error: ' num2str(min_err(1))]; ...
                    ['Curve Number: ' num2str(c_n)]} ;
            text(0.1,0.2,information);

            set(gcf,'Position',[30 60 880 500]);
        end

        total_error{c_n,1} = total_err;
        curves{c_n,5} = FC;
        curves{c_n,6} = E_pointwise;
        curves{c_n,7} = E_apparent;
        clear total_err E_apparent E_pointwise FC
        clear extension retract pre_* post_* min_err
    end
end

end
