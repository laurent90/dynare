function a=set_Kalman_starting_values(a,M_,oo_,options_,bayestopt_)
% function a=set_Kalman_starting_values(a,M_,oo_,options_,bayestopt_)
% Sets initial states guess for Kalman filter/smoother based on histval 
% 
% INPUTS 
%   o a             [double]   (p*1) vector of states
%   o M_            [structure] decribing the model
%   o oo_           [structure] storing the results
%   o options_      [structure] describing the options
%   o bayestopt_    [structure] describing the priors
%  
% OUTPUTS
%   o a             [double]    (p*1) vector of set initial states

if ~isempty(M_.endo_histval)
    if options_.loglinear && ~options_.logged_steady_state
        a(bayestopt_.mf0)      = log(M_.endo_histval(oo_.dr.order_var(oo_.dr.restrict_var_list(bayestopt_.mf0)),:)) - log(oo_.dr.ys(oo_.dr.order_var(oo_.dr.restrict_var_list(bayestopt_.mf0))));            
    elseif ~options_.loglinear && ~options_.logged_steady_state
        a(bayestopt_.mf0)      = M_.endo_histval(oo_.dr.order_var(oo_.dr.restrict_var_list(bayestopt_.mf0)),:) - oo_.dr.ys(oo_.dr.order_var(oo_.dr.restrict_var_list(bayestopt_.mf0)));        
    else
        error(['The steady state is logged. This should not happen. Please contact the developers'])
    end
end
