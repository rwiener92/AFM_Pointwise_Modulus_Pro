function E_INT_pointwise = E_DEPTH_LININTERP(E_pointwise,min_depth,max_depth,stp_depth)
% Designed for call from afm_post_pro_editor
% inputs depth dependent moduli data calculated from @calculate_moduli
% Returns E_pointwise matrix that has been discretized into symetrically
% spaced data via linear interpolation for comparison across multiple
% data sets
%
% ------------------------------------------------------------------------

depth = E_pointwise(:,1); % nanometers
mod_depth = E_pointwise(:,2); % kilopascals

m = 1;
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
            E_INT_pointwise(e_count,2) = mod_depth(m) + (k - depth(m))*...
                ((mod_depth(m+1)-mod_depth(m))/(depth(m+1)-depth(m)));
        end     
    else
        E_INT_pointwise(e_count,2) = NaN;
        E_INT_pointwise(e_count,1) = NaN;
    end
E_INT_pointwise(e_count,1) = k;

end

end