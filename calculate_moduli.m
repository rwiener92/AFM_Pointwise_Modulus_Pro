function [E_apparent,E_pointwise] = calculate_moduli(...
    FC,geometry_opt,pointwise_opt,display_opt)
% CALCULATE_MODULI calculates the depth dependent and average apparent
%   moduli using the supplied force-depth curve and geometry information.
%   E_APPARENT = CALCULATE_MODULI(FC,GEOMETRY_OPT) uses the geometry
%       assumption defined by the scalar GEOMETRY_OPT and calculates the
%       apparent modulus by fitting the appropriate function to the
%       force-depth data.
%   [E_APPARENT,E_POINTWISE] =
%       CALCULATE_MODULI(FC,GEOMETRY_OPT,POINTWISE_OPT) also calculates the
%       pointwise modulus as a function of depth with the resolution
%       defined by the scalar POINTWISE_OPT.
%   [E_APPARENT,E_POINTWISE] =
%       CALCULATE_MODULI(FC,GEOMETRY_OPT,POINTWISE_OPT,DISPLAY_OPT)
%       displays the computational progress depending on the value of the
%       DISPLAY_OPT.
%
%   GEOMETRY_OPT - 1: Standard Si3N4 tip
%                  2: Spherical tip with radius R
%                  3: Cylindrical (flat) tip with radius R
%                  4: Custom tip
%   POINTWISE_OPT- 0: Do not calculate depthwise modulus
%                  1: Calculate depthwise modulus at each depth point
%                  2: Calculate depthwise modulus at 50nm increments
%   DISPLAY_OPT -  0: Do not display progress
%                  1: Display depthwise calculation count (n/n_total)
%                  2: Display the current depth
%
%   AZELOGLU E.U. (c) 2005

% Note: Postmark option on line 44 is removed; therefore 50nm increments do
% not work at the time.
if nargin == 3
    display_opt = 0;
elseif nargin == 2
    display_opt = 0;
    pointwise_opt = 0;
end

if nargout < 2
    pointwise_opt = 0;
end

post_sample_no = size(FC,1);
% if pointwise_opt == 1
    post_mark = [1:1:post_sample_no];
% else
%     for n = 2:post_sample_no
%         % Below is a quick code that marks the locations (in the
%         % force-depth graph) that are near multiples of fifty. The code
%         % will remember these locations and if it's asked to operate faster
%         % (i.e. calculate pointwise moduli at discretized depths instead
%         % every single number) then it would look at the depth markers that
%         % it has stored and go to those points to evaluate moduli much
%         % faster. EXAMPLE: If the depth is higher than fifty at the point
%         % of interest and lower than fifty at the previous location, then
%         % it would store both of the points (since chances of a perfect
%         % match is nearly none, that is not considered in this code). Below
%         % each odd column of the post_mark vector is a value that marks the
%         % 50 increment, the following even column is the followup value
%         % (i.e. 49 and 51):
%         k = 1;
%         for m = 50:50:800
%             if (FC(n,1) > m) & (FC(n-1,1) < m)
%                 post_mark(k) = n-1;
%                 post_mark(k+1) = n;
%             end
%             k = k + 2;
%         end
% 	end
% end
% During the pointwise modulus (E_pointwise) calculation the parameters at
% the right hand side of the force-depth equation are clustered as the
% x-variable (x_for_E_fit) and force is clustered as the y-variable
% (y_for_E_fit), which will be used to calculate an apparent modulus
% (E_apparent) for the given contact geometry by fitting a least squares to
% y = E*x equation.


% -------------------------------------------------------------------------
% MAIN MODULUS CALCULATION KERNEL
% -------------------------------------------------------------------------
warning off MATLAB:divideByZero;
if geometry_opt == 1
    % OPTION 1 is the standard blunt cone Si3N4 tip with previously
    % supplied cantilever stiffness.
	R_tip = 42; % Radius of curvature
	t_tip = 40*pi/180; % Cone semi-angle
	b_tip = R_tip*t_tip; % Cylindrical radius
    
    k = 1;
    for n = post_mark
        if isnan(FC(n,1)) ~= 1
            h_tip = FC(n,1);
            a_tip = sqrt(R_tip * h_tip);
            % The solution for the modulus uses a spherical geometry until
            % the contact radius surpasses the cylindrical radius (b_tip)
            % and a blunt-cone geometry once this point is reached:
            if a_tip < b_tip
                if pointwise_opt ~= 0
                    E_temp(k,1) = h_tip;
                    E_temp(k,2) = (0.75*FC(n,2))/(sqrt(R_tip)*h_tip^1.5);
                end
                x_for_E_temp(k) = (4/3)*sqrt(R_tip)*h_tip^1.5;
            elseif a_tip >= b_tip
                % When the critical depth is reached, contact radius will
                % surpass the spherical blunt tip radius. At that point,
                % summon the GET_CONTACT_RADIUS function that numerically
                % solves for the contact radius using equation (4) of
                % Briscoe et al.
                a_tip = get_contact_radius(b_tip,h_tip,R_tip,t_tip);
                messy_part = (a_tip*h_tip - (a_tip^2/(2*tan(t_tip)))*...
                    (0.5*pi - asin(b_tip/a_tip)) - (a_tip^3/(3*R_tip))...
                    + sqrt(a_tip^2 - b_tip^2)*((b_tip/(2*tan(t_tip)))...
                    + (a_tip^2 - b_tip^2)/(3*R_tip)));
                if pointwise_opt ~= 0
                    E_temp(k,1) = h_tip;
                    E_temp(k,2) = (0.5*FC(n,2))/messy_part;
                end
                x_for_E_temp(k) = 2*messy_part;
            end
            y_for_E_temp(k) = FC(n,2);
            a_vector(k) = a_tip;
            h_vector(k) = h_tip;
            k = k + 1;
        end
    end
    
    if pointwise_opt == 1
        E_pointwise = E_temp;
        x_for_E_fit = x_for_E_temp;
        y_for_E_fit = y_for_E_temp;
        clear *temp
    else
        k = 1;
        for n = 1:(0.5*length(post_mark))
            depth_gap = FC(post_mark(k+1),1) - FC(post_mark(k),1);
            interp_constant = n*50 - FC(post_mark(k),1);
            x_for_E_fit(n) = x_for_E_temp(k) + (interp_constant * ...
                ((x_for_E_temp(k+1) - x_for_E_temp(k))/depth_gap));          
            y_for_E_fit(n) = y_for_E_temp(k) + (interp_constant * ...
                ((y_for_E_temp(k+1) - y_for_E_temp(k))/depth_gap));        
            if pointwise_opt == 2
                E_pointwise(n,1) = E_temp(k,1) + (interp_constant *...
                    ((E_temp(k+1,1) - E_temp(k,1))/depth_gap));  
                E_pointwise(n,2) = E_temp(k,2) + (interp_constant *...
                    ((E_temp(k+1,2) - E_temp(k,2))/depth_gap));
            end
            k = k + 2;
        end
    end  
    E_apparent = sum(y_for_E_fit .* x_for_E_fit)/sum(x_for_E_fit .^ 2);  
    
    
    
    
% -------------------------------------------------------------------------
elseif geometry_opt == 2
    % OPTION 2 is the standard spherical probe with user supplied radius
    % and previously supplied cantilever stiffness.
    R_tip = 7500; % Make this a variable
    x_for_E_fit = (4/3)*sqrt(R_tip)*(FC(:,1).^1.5);
    y_for_E_fit = FC(:,2);
    E_apparent = sum(y_for_E_fit .* x_for_E_fit)/sum(x_for_E_fit .^ 2);
    if pointwise_opt ~= 0
        for n = post_mark
            E_pointwise(n,1) = FC(n,1);
            E_pointwise(n,2) = (0.75*FC(n,2))/(sqrt(R_tip)*FC(n,1)^1.5);
        end
    end 
    
    
    
    
% -------------------------------------------------------------------------    
elseif geometry_opt == 3
    % OPTION 3 is the standard cylindrical (hemispherical) probe with user
    % supplied radius and previously supplied cantilever stiffness.
    b_tip = 7500; % Make this a variable
    x_for_E_fit = 2*b_tip*FC(:,1);
    y_for_E_fit = FC(:,2);
    E_apparent = sum(y_for_E_fit .* x_for_E_fit)/sum(x_for_E_fit .^ 2);
    if pointwise_opt ~= 0
        for n = post_mark
            if pointwise_opt == 1
                E_pointwise(n,1) = FC(n,1);
                E_pointwise(n,2) = (0.5*FC(n,2))/(b_tip * FC(n,1));
            elseif pointwise_opt == 2
                E_temp(n,1) = FC(n,1);
                E_temp(n+1,1) = FC(n+1,1);
            end
            
            if display_opt == 1
                disp(['Computing pointwise modulus... (' num2str(n) ...
                        '/' num2str(post_sample_no) ')']);
            elseif display_opt == 2
                disp(['Computing modulus at depth: ' num2str(FC(n,1)) '']);
            end
        end
    end
    
    
    
% -------------------------------------------------------------------------    
elseif geometry_opt == 4
    % OPTION 4 is the custom operation mode which asks for user input for
    % each parameter.
	R_tip = 50; % Make this a variable
	t_tip = 40*pi/180; % Make this a variable
    smooth_transition = 'Y'; % Make this a variable

    if smooth_transition == 'Y' | smooth_transition == 'y'
        b_tip = R_tip*t_tip;
    else
        b_tip = 30; % Make this a variable
    end
    
    for n = 1:post_sample_no
        h_tip = FC(n,1);
        % Summon the GET_CONTACT_RADIUS function that numerically solves
        % for the contact radius using equation (4) of Briscoe et al.
        a_tip = get_contact_radius(b_tip,h_tip,R_tip,t_tip);
        % The solution for the modulus uses a spherical geometry until the
        % contact radius surpasses the cylindrical radius (b_tip) and a
        % blunt-cone geometry once this point is reached:
        if a_tip <= b_tip
            a_tip = sqrt(R_tip * h_tip);
            E_pointwise(n,1) = h_tip;
            E_pointwise(n,2) = (0.75*FC(n,2))/(sqrt(R_tip)*h_tip^1.5);
            x_for_E_fit(n) = (4/3)*sqrt(R_tip)*h_tip^1.5;
        elseif a_tip > b_tip
            messy_part = (a_tip*h_tip - (a_tip^2/(2*tan(t_tip)))*(0.5*pi - ...
                asin(b_tip/a_tip)) - (a_tip^3/(3*R_tip)) + sqrt(a_tip^2 - ...
                b_tip^2)*((b_tip/(2*tan(t_tip))) + (a_tip^2 - b_tip^2)/...
                (3*R_tip)));
            E_pointwise(n,1) = h_tip;
            E_pointwise(n,2) = (0.5*FC(n,2))/messy_part;
            x_for_E_fit(n) = 2*messy_part;
        end
        y_for_E_fit(n) = FC(n,2);
        a_vector(n) = a_tip;
        h_vector(n) = h_tip;     
    end  
    E_apparent = sum(y_for_E_fit .* x_for_E_fit)/sum(x_for_E_fit .^ 2);  
    
    
end


% Convert the moduli values to kilopascals (kPa):
E_pointwise(:,2) = E_pointwise(:,2) * 1e6;
E_apparent = E_apparent * 1e6;
warning on all;
