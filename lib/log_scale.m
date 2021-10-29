function [x_log] = log_scale(x,t,tr,st)
%log_scale(variable, time array, dt, stopping time)
%   This function takes an input of two arrays, the first being the
%   variable of interest and the second being the time array for that
%   variable. Using a while loop, it resamples the arrays into logarithmic
%   form. The final setting chooses whether or not to skip negative values
    
    n_log = length(x);
    n_logMax = 10;
    i = 1;
    j = 1;
    x_log = [];
    tr_loop = tr;
    
    while i < n_log
        % How many points do we have in this decade?
        magPoints = (find( floor(log10(t)) == floor(log10(j*tr_loop)) ));
        if length(magPoints) < 10
            stepsize = 1;
        else
            stepsize = median( floor( diff( linspace(magPoints(1), magPoints(end), 10) ) ) );
        end

        % If you are between the first and last time...
        if ((t(i) >= j*tr_loop) || (abs(t(i) - j*tr_loop) < eps)) && t(i) <= st
            % Note in this logical check that the first true/false
            % evaluation has to contain logica that checks if we are
            % basically "in the same bin" (i.e. the second check for the
            % difference between numbers being less than the machine
            % error). If you don't check this difference, you can
            % erroneously skip datapoints that should, indeed, be included
            % because of machine error.
            
            x_log = [x_log x(i)];
            j = j + 1;
            
            % If you've hit the top of this log range...
            if j == n_logMax
                tr_loop = tr_loop * 10;
                j = 1;
            end
        end

        i = i + stepsize;
        
    end

end

