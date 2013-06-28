function [SE] = get_Standard_Errors_SMM(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,SMMinfo,DynareResults)
% [SE] = get_Standard_Errors_SMM(xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,SMMinfo,DynareResults)
% This function computes standard errors to the SMM estimates
% By Martin Andreasen

% INPUTS 
%   o xparam1:                  initial value of estimated parameters as returned by set_prior()
%   o DynareDataset:            data after required transformation
%   o DynareOptions             Matlab's structure describing the options (initialized by dynare, see @ref{options_}).
%   o Model                     Matlab's structure describing the Model (initialized by dynare, see @ref{M_}).          
%   o EstimatedParameters:      Matlab's structure describing the estimated_parameters (initialized by dynare, see @ref{estim_params_}).
%   o SMMInfo                   Matlab's structure describing the SMM settings (initialized by dynare, see @ref{bayesopt_}).
%   o DynareResults             Matlab's structure gathering the results (initialized by dynare, see @ref{oo_}).
%  
% OUTPUTS 
%   o SE                       [nparam x 1] vector of standard errors
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 2013 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.


% Get the Jacobian of the moment difference function
D = fdjac('get_SMM_objectFun_moments',xparam1,DynareDataset,DynareOptions,Model,EstimatedParameters,SMMinfo,DynareResults);

W=DynareResults.smm.W;
Wopt=DynareResults.smm.Wopt;

T    = DynareDataset.info.ntobs; %Number of observations

if ~DynareOptions.smm.optimal_weighting
    AVar = 1/T*DynareResults.smm.variance_correction_factor*eye(length(xparam1),length(xparam1))*((D'*W*D)\D'*W/Wopt*W*D/(D'*W*D));
else
    AVar = 1/T*DynareResults.smm.variance_correction_factor*eye(length(xparam1),length(xparam1))/(D'*DynareResults.smm.W*D);
end
SE   = sqrt(diag(AVar));