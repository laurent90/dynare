function [scale,shape] = weibull_specification(mu,sigma) % --*-- Unitary tests --*--
% Computes the Weibull hyperparameters from the prior mean and standard deviation.
%
%@info:
%! @deftypefn {Function File} {[@var{scale}, @var{shape} ]=} colon (@var{mu}, @var{sigma}
%! @anchor{distributions/weibull_specification}
%! @sp 1
%! Computes the Weibull hyperparameters from the prior mean (@var{mu}) and standard deviation (@var{sigma}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item mu
%! Double scalar, prior mean.
%! @item sigma
%! Positive double scalar, prior standard deviation.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item scale
%! Positive double scalar, first hypermarameter (scale) of the Weibull prior .
%! @item shape
%! Positive double scalar, second hypermarameter (shape) of the Weibull prior.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{set_prior}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! 
%! @end deftypefn
%@eod:

% Copyright (C) 2015 Dynare Team
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
% You should have received scale copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

check_solution_flag = 1;
scale = [];
shape = [];

if sigma^2 < Inf
    options=optimset('MaxFunEvals',10000,'TolF',1e-10,'Display','off');
    x = fsolve(@(x)weibull_specification_error(x,mu,sigma),[mu;sigma],options);
    scale=exp(x(1));
    shape=exp(x(2));
    if check_solution_flag
        if abs(mu-scale*gamma(1+1/shape))>1e-7
            error('weibull_specification:: Failed in solving for the hyperparameters!');
        end
        if abs(sigma^2-scale^2*(gamma(1+2/shape)-(gamma(1+1/shape))^2))>1e-7
            error('weibull_specification:: Failed in solving for the hyperparameters!');
        end
    end
else
    error('weibull_specification:: Infinite variance not allowed for Weibull prior')
end

end

function outvalue=weibull_specification_error(x,mu,sigma)
    scale_par=exp(x(1,:)); %make sure only positive solutions are returned
    shape_par=exp(x(2,:)); %make sure only positive solutions are returned
    mean_gamma=scale_par.*gamma(1+1./shape_par);
    variance_gamma=scale_par.^2.*(gamma(1+2./shape_par)-(gamma(1+1./shape_par)).^2);
    outvalue(1,:)=mean_gamma-mu;
    outvalue(2,:)=variance_gamma-sigma.^2;
end

%@test:1
%$
%$ [a1,b1] = weibull_specification(1,1)
%$ [a2,b2] = weibull_specification(0.01,5)
%$ % Check the results.
%$ t(1) = dassert(a1,1,1e-6);
%$ t(2) = dassert(b1,1,1e-6);
%$ t(3) = dassert(a2,1.6156e-09,1e-6);;
%$ t(4) = dassert(b2,0.0978,1e-6);
%$ T = all(t);
%@eof:1