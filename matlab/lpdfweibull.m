function  [ldens,Dldens,D2ldens] = lpdfweibull(x,a,b)
% Evaluates the logged UNIVARIATE WEIBULL PDF at x.
%
% INPUTS 
%    x     [double]  m*n matrix of locations,
%    a     [double]  m*n matrix or scalar, First WEIBULL distribution parameters (scale) 
%    b     [double]  m*n matrix or scalar, Second WEIBULL distribution parameters (shape). 
%
% OUTPUTS 
%    ldens [double]  m*n matrix of logged WEIBULL densities evaluated at x.
%     
%        
% SPECIAL REQUIREMENTS
%    none

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
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

ldens = log(wblpdf(x,a,b));

ldens = -Inf(size(x)) ;
idx = find(x>0);

if length(a)==1
    ldens(idx) = log(b)-log(a)+(b-1).*(log(x(idx))-log(a))-(x(idx)./a).^b; 
else
    ldens(idx) = log(b(idx))-log(a(idx))+(b(idx)-1).*(log(x(idx))-log(a(idx)))-(x(idx)./a(idx)).^b(idx); 
end

if nargout >1 
    if length(a)==1
        Dldens(idx) =  (b-1 - b.*(x(idx)./a).^b)./x(idx); %(d)/(dx)(log((b (x/a)^(b-1) exp(-(x/a)^b))/a)) = (-b (x/a)^b+b-1)/x
    else
        Dldens(idx) =  (b(idx)-1 - b(idx).*(x(idx)./a(idx)).^b(idx))./x(idx); %(d)/(dx)(log((b (x/a)^(b-1) exp(-(x/a)^b))/a)) = (-b (x/a)^b+b-1)/x
    end
end

if nargout == 3 
    if length(a)==1
        D2ldens(idx) =  -((b-1).*(b.*(x(idx)./a).^b+1))./x(idx).^2; %(d)/(dx)((-b (x/a)^b+b-1)/x) = -((b-1) (b (x/a)^b+1))/x^2
    else
        D2ldens(idx) =  -((b(idx)-1).*(b(idx).*(x(idx)./a(idx)).^b(idx)+1))./x(idx).^2; %(d)/(dx)((-b (x/a)^b+b-1)/x) = -((b-1) (b (x/a)^b+1))/x^2
    end
end