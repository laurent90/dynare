function w_rnd = wblrnd(scale,shape,varargin) % --*-- Unitary tests --*--
% function w_rnd = wblrnd(scale,shape,varargin)
% Return a matrix of random samples from the Weibull distribution with
% parameters @var{scale} and @var{shape}.
%
% When called with a single size argument, return a square matrix with
% the dimension specified.  When called with more than one scalar argument the
% first two arguments are taken as the number of rows and columns and any
% further arguments specify additional matrix dimensions.  The size may also
% be specified with a vector of dimensions @var{sz}.
%
% If no size arguments are given then the result matrix is the common size of
% @var{scale} and @var{shape}.


% Adapted for Matlab (R) from GNU Octave 4.0
% Original file: statistics/distributions/wblrnd.m
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

if nargin<2
    error(message('wblrnd: Note enough input arguments'));
end

if ~isequal(size(scale),size(shape))
    error ('wblrnd: SCALE, and SHAPE must be of common size or scalars');
end

if (any(~isreal(scale)) || any(~isreal(shape)))
    error ('wblrnd: SCALE and SHAPE must not be complex');
end

if (nargin == 2)
    sz = size (scale);
elseif (nargin == 3)
    if (isscalar(varargin{1}) && varargin{1} >= 0)
        sz = [varargin{1}, varargin{1}];
    elseif (isrow (varargin{1}) && all (varargin{1} >= 0))
        sz = varargin{1};
    else
        error ('wblrnd: dimension vector must be row vector of non-negative integers');
    end
elseif (nargin > 3)
    if (any (cellfun (@(x) (~isscalar(x) || x < 0), varargin)))
        error ('wblrnd: dimensions must be non-negative integers');
    end
    sz = [varargin{:}];
end

if (~isscalar (scale) && ~isequal(size(scale),sz))
    error ('wblrnd: SCALE and SHAPE must be scalar or of size SZ');
end


if (isscalar (scale) && isscalar (shape))
    if ((scale > 0) && (scale < Inf) && (shape > 0) && (shape < Inf))
        w_rnd = scale*(-log(rand(sz))).^(1/shape);
    else
        w_rnd = NaN (sz);
    end
else
    w_rnd = scale.*(-log(rand(sz))).^(1./shape);
    
    k = (scale <= 0) | (scale == Inf) | (shape <= 0) | (shape == Inf);
    w_rnd(k) = NaN;
end

end

%@test:1
%$
%$ t(1)=dassert (size (wblrnd (1,2)), [1, 1]);
%$ t(2)=dassert (size (wblrnd (1, 2, 3)), [3, 3]);
%$ t(3)=dassert (size (wblrnd (1, 2, [4 1])), [4, 1]);
%$ t(4)=dassert (size (wblrnd (1, 2, 4, 1)), [4, 1]);
%$ r_vec=wblrnd(1, 2, 1000000, 1); %Wolfram Alpha: WeibullDistribution[2,1]
%$ t(5)=dassert(mean(r_vec),sqrt(pi)/2,1e-3);
%$ t(6)=dassert(std(r_vec),sqrt(1-pi/4),1e-3);
%$ T = all(t);
%$
%@eof:1