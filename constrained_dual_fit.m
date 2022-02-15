function [a_coefficients,total_rms,convergence_flag] = ...
    constrained_dual_fit(pre_fit,post_fit,pre_n,post_n)
% CONSTRAINED_DUAL_FIT calculates the total error for a two-phase discrete
%   polynomial fit with an arbitrary contact point prescribed. All the
%   points before the prescribed point belongs to the pre-contact curve,
%   whereas all points after belong to the post-contact portion. The
%   post-contact function uses x-x* as the independent variable, where x*
%   is the prescribed contact coordinate.
%
%   A_COEFFICIENTS = CONSTRAINED_DUAL_FIT(PRE_FIT,POST_FIT) computes the
%   optimum polynomial coefficients, A_COEFFICIENTS for pre-contact
%   (PRE_FIT) and post-contact (POST_FIT)curves using a quadratic-cubic
%   fits, respectively. It constrains the contact point y-value for both
%   functions to be equal as in: y* = a + bx* + cx*^2 = c. PRE_FIT and
%   POST_FIT are column matrices with the first column being the
%   independent variable and the second one its dependent.
%
%   A_COEFFICIENTS = CONSTRAINED_DUAL_FIT(PRE_FIT,POST_FIT,PRE_N,POST_N)
%   introduces two optional power terms that limit the polynomials orders
%   to smaller values, such as linear-quadratic or linear-linear. AFM Post
%   Processor code by default uses linear-quadratic.
%
%   AZELOGLU E.U. (c) 2005
%   Main kernel adapted from COSTA K.D.

if nargin == 2
    pre_n = 2;
    post_n = 3;
elseif nargin == 3
    post_n = 3;
end

if pre_n > 2 | post_n > 3
    error('Cannot fit higher than quadratic-cubic polynomials');
end

% -------------------------------------------------------------------------
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
xs = post_fit(1,1); % Prescribed contact point, x_star

n1 = size(pre_fit,1);   % Pre-contact number of variables
x1 = pre_fit(:,1);      % Pre-contact x-vector
y1 = pre_fit(:,2);      % Pre-contact y-vector

n2 = size(post_fit,1);  % Post-contact number of variables
x2 = (post_fit(:,1)-xs);     % Post-contact x-vector
% Convert x2 into x - x* so that the constraint is easier to apply.
y2 = post_fit(:,2);     % Post-contact y-vector

options = optimset('display','off','LargeScale','on',...
    'Algorithm','active-set','TolX',1e-10);
xdata = [x1',x2'];
ydata = [y1',y2'];

% First compute the unconstrained functions' coefficients:
pre_temp = polyfit(x1,y1,pre_n);
post_temp = polyfit(x2,y2,post_n);
% Then plug them in as the initial guesses for the constrained functions'
% coefficients:
guess = [pre_temp,post_temp];
% Herein the first "PRE_N + 1" parameters belong to the pre-contact function
% whereas the latter "POST_N + 1" parameters belong to the post-contact
% curve.

% A secondary way to apply the equality constraint is through the Aeq and
% beq arguments that might be introduces to the constraint module:
% Aeq = [1 pre_fit(:,1) pre_fit(:,1).^2];
% beq = guess(4);

a_coefficients = fmincon(@dual_poly,guess,[],[],[],[],[],[],...
    @equality_con,options,xdata,ydata,n1,pre_n,post_n);
if nargout == 2
	[a_coefficients total_rms] = fmincon(@dual_poly,guess,[],[],[],[],...
        [],[],@equality_con,options,xdata,ydata,n1,pre_n,post_n);
elseif nargout == 3
	[a_coefficients total_rms convergence_flag] = fmincon(@dual_poly,guess,...
        [],[],[],[],[],[],@equality_con,options,xdata,ydata,n1,pre_n,post_n);
end    

% !! For Debugging Purposes Only !!
% Evaluate and plot the fits for visual inspection:
% pre_poly = polyval(a_coefficients(1:pre_n+1),x1);
% post_poly = polyval(a_coefficients(pre_n+2:end),x2);
% plot(x1,y1,'.b',x1,pre_poly,'-b',x2+xs,y2,'.r',x2+xs,post_poly,'-r');


function F = dual_poly(guess,xdata,ydata,n1,pre_n,post_n)
for n = 1:length(xdata)
    if n <= n1
        if pre_n == 2
            f(n) = (ydata(n) - (guess(1)*xdata(n)^2 + guess(2)*xdata(n) + ...
                guess(3)))^2;
        elseif pre_n == 1
            f(n) = (ydata(n) - (guess(1)*xdata(n) + guess(2)))^2;
        end
    elseif n > n1
        if post_n == 3
            f(n) = (ydata(n) - (guess(pre_n+2)*xdata(n)^3 + ...
                guess(pre_n+3)*xdata(n)^2 + guess(pre_n+4)*xdata(n) ...
                + guess(pre_n+5)))^2;
        elseif post_n == 2
            f(n) = (ydata(n) - (guess(pre_n+2)*xdata(n)^2 + ...
                guess(pre_n+3)*xdata(n) + guess(pre_n+4)))^2;
        elseif post_n == 1
            f(n) = (ydata(n) - (guess(pre_n+2)*xdata(n) + guess(pre_n+3)))^2;
        end
    end
end
F = sum(f(:));
warning('on','all');

function [c,ceq] = equality_con(a_coefs,xdata,ydata,n1,pre_n,post_n)
guess = a_coefs;
c = [];
if pre_n == 2
    ceq = guess(1)*xdata(n1)^2 + guess(2)*xdata(n1) + ...
        guess(3) - guess(end);
elseif pre_n == 1
    ceq = guess(1)*xdata(n1) + guess(2) - guess(end);
end

