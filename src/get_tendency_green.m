function [dCdt] = get_tendency_green(t,C,L)
% function [dCdt] = get_tendency_green(t,C,L)

    verbose = true;
    % does is slow down the code to execute an `if` statement every timestep?
    if verbose
        display(['get tendency: t=',num2str(t),' yr'])
    end
    %t
    dCdt = L*C; 
end
