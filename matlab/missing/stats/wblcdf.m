function w_cdf = wblcdf(x,scale,shape) % --*-- Unitary tests --*--
% function pdf = wblcdf(x,scale,shape)
% Compute the cumulative distribution function (CDF) at x of the Weibull distribution with scale parameter scale and shape parameter shape.
%
% This is defined as
%   1 - exp (-(x/scale)^shape)
% for x >= 0.
% 
% Default values are @var{scale} = 1, @var{shape} = 1.


% Adapted for Matlab (R) from GNU Octave 4.0
% Original file: statistics/distributions/wblcdf.m
% Original author: KH <Kurt.Hornik@wu-wien.ac.at>

% Copyright (C) 2012 Rik Wehbring
% Copyright (C) 1995-2015 Kurt Hornik
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

if nargin<1
    error(message('wblcdf: Note enough input arguments'));
end
if nargin < 2
    scale = ones(size(x));
end
if nargin < 3
    shape = ones(size(x));
end

if ~(isscalar(scale) && isscalar(shape)) && ~isequal(size(scale),size(shape),size(x))
    error ('wblcdf: X, SCALE, and SHAPE must be of common size or scalars');
end

if (any(~isreal(x)) || any(~isreal(scale)) || any(~isreal(shape)))
    error ('wblcdf: X, SCALE, and SHAPE must not be complex');
end

w_cdf = NaN(size(x));

ok_condition = ((scale > 0) & (scale < Inf) & (shape > 0) & (shape < Inf));

k = (x <= 0) & ok_condition;
w_cdf(k) = 0;

k = (x == Inf) & ok_condition;
w_cdf(k) = 1;

k = (x > 0) & (x < Inf) & ok_condition;
if (isscalar(scale) && isscalar(shape))
    w_cdf(k) = 1 - exp(-(x(k)/scale).^shape);
  else
    w_cdf(k) = 1 - exp(-(x(k)./scale(k)).^shape(k));
end

end

%@test:1
%$
%$ x = [-1 0 0.5 1 Inf]';
%$ y = [0; 1-exp(-x(2:4)); 1];
%$ t(1) = dassert(wblcdf (x,ones(5,1),ones (5,1)),y);
%$ t(2) = dassert(wblcdf (x,1,1),y);
%$ T = all(t);
%$
%@eof:1