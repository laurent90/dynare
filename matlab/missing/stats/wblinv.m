function invcdf = wblinv(x,scale,shape) % --*-- Unitary tests --*--
% function pdf = wblinv(x,scale,shape)
% Compute the quantile (the inverse of the CDF) at @var{x} of the
% Weibull distribution with scale parameter @var{scale} and
% shape parameter @var{shape}.
%
% Default values are @var{scale} = 1, @var{shape} = 1.


% Adapted for Matlab (R) from GNU Octave 4.0
% Original file: statistics/distributions/wblinv.m
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
    error(message('wblinv: Note enough input arguments'));
end
if nargin < 2
    scale = ones(size(x));
end
if nargin < 3
    shape = ones(size(x));
end

if ~(isscalar(scale) && isscalar(shape)) && ~isequal(size(scale),size(shape),size(x))
    error ('wblinv: X, SCALE, and SHAPE must be of common size or scalars');
end

if (any(~isreal(x)) || any(~isreal(scale)) || any(~isreal(shape)))
    error ('wblinv: X, SCALE, and SHAPE must not be complex');
end

invcdf = NaN (size (x));

ok_condition = ((scale > 0) & (scale < Inf) & (shape > 0) & (shape < Inf));

k = (x == 0) & ok_condition;
invcdf(k) = 0;

k = (x == 1) & ok_condition;
invcdf(k) = Inf;

k = (x > 0) & (x < 1) & ok_condition;
if (isscalar (scale) && isscalar (shape))
    invcdf(k) = scale*(-log(1-x(k))).^(1/shape);
else
    invcdf(k) = scale(k).*(-log(1-x(k))).^(1./shape(k));
end

%@test:1
%$
%$ x = [-1 0 0.63212055882855778 1 2];
%$ t(1)=dassert(wblinv(x,ones(1,5),ones(1,5)),[NaN 0 1 Inf NaN],eps);
%$ t(2)=dassert(wblinv(x,1,1),[NaN 0 1 Inf NaN],eps);
%$ T = all(t);
%$
%@eof:1