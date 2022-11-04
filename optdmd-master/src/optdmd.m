function [w,e,b,varargout] = optdmd(X,t,r,imode,varargin)
%OPTDMD Wrapper of VARPRO2 for computing the optimized DMD of data
%
%   [w,e,b,varargout] = optdmd(X,t,r,imode,varargin)
%
% Input:
%
% X - data matrix, length(t) columns with each a snapshot
%   of some time-varying system
% t - times corresponding to each snapshot
% r - rank of fit, i.e. number of exponentials to fit
% imode - flag, determines type of computation
%   imode = 1, fit full data, slower
%   imode = 2, fit data projected onto first r POD modes
%      or columns of varargin{2} (should be at least r
%      columns in varargin{2})
% varargin{1} - options structure for varpro2. see varpro_opts.m
%   for details
% varargin{2} - initial guess. if not provided, one is computed
%   using trapezoidal rule approximation
% varargin{3} - orthogonal basis for projection (if POD modes precomputed
%   or a different basis desired)
% varargin{4} - linear constraint options structure.
%                       See varpro_lsqlinopts.m for details.
%                       allows you to enforce linear constraints.
% varargin{5} - gamma = tikhonov regularization term. if provided
%                       the minimization problem becomes 
%
%                min  | y - phi*b |_F^2 + | gamma alpha |_2^2 
%
%               where gamma is either a scalar or matrix.  
%
% varargin{6} = prox - proximal operator to be applied to the 
%                      vector alpha at each step (e.g. projection
%                      onto a set). For now, this overrides the 
%                      lsqlin options in varargin{4}.
%
%
% Output:
%
% w - each column is a DMD mode
% e - each entry e(i) is an eigenvalue corresponding to w(:,i)
% b - the best fit coefficient of each DMD mode
% varargout{1} - return projected system matrix
%   A =  (u'*w)*diag(e)*(pinv(w)*u) where u is the first
%   r POD modes or varargin{3}.
% varargout{2} - return basis for projection
% varargout{3} - return full system matrix A = w*diag(e)*pinv(w)
%
% X should be approximated by
%
% X ~ w*diag(b)*exp(e*t')
%
% if t is a column vector.
%
% Examples:
%
%   >> [w,e,b] = optdmd(X,t,r,imode);
%   >> [w,e,b] = optdmd(xdata,ts,r,imode,opts,[],u);
%   >> [w,e,b] = optdmd(xdata,ts,r,imode,[],e_init);
%   >> [w,e,b,atilde,u,afull] = optdmd(X,t,r,imode);
%   >> [w,e,b] = optdmd(xdata,ts,r,imode,[],e_init,[],copt);
%
% See also VARPRO_OPTS, VARPRO2

%
% Copyright Travis Askham 2017
%
% MIT License
%

if (nargin >= 7 && ~isempty(varargin{3}))
    u = varargin{3};
elseif ( imode==2 || nargout > 3 || nargin < 6 || isempty(varargin{2}) )
    [u,~,~] = svd(X,'econ');
    if (imode == 2)
        u = u(:,1:r);
    end
end

if (nargin < 5 || isempty(varargin{1}))
    opts = varpro_opts();
else
    opts = varargin{1};
end

if (nargin < 6 || isempty(varargin{2}))
    
    % use projected trapezoidal rule approximation
    % to eigenvalues as initial guess
    
    ux1 = u'*X;
    ux2 = ux1(:,2:end);
    ux1 = ux1(:,1:end-1);
    
    t1 = t(1:end-1);
    t2 = t(2:end);
    
    dx = (ux2-ux1)*diag(1./(t2-t1));
    xin = (ux1+ux2)/2;
    
    [u1,s1,v1] = svd(xin,'econ');
    
    u1 = u1(:,1:r);
    v1 = v1(:,1:r);
    s1 = s1(1:r,1:r);
    
    atilde = u1'*dx*v1/s1;
    
    alpha_init = eig(atilde);
    
    clear ux1 ux2 atilde t1 t2 dx xin
    
else
    
    % use user provided initial guess
    
    alpha_init = varargin{2};
end

% check if linear constraints are supplied

if (nargin < 8 || isempty(varargin{4}))
    copts = [];
else
    copts = varargin{4};
end

% check if Tikhinov regularization is supplied

if (nargin < 9 || isempty(varargin{5}))
    gamma = [];
else
    gamma = varargin{5};
end

% check if proxfun is supplied 

if (nargin < 10 || isempty(varargin{6}))
    proxfun = [];
else
    proxfun = varargin{6};
end

if (imode == 2)
    
    % projected version
    
    m = length(t);
    [~,n] = size(u);
    ia = r;
    is = r;
    [w,e,~,~,~,~] = varpro2(transpose(u'*X),t, ...
            @varpro2expfun,@varpro2dexpfun,m,n,is,ia,alpha_init, ...
            opts,copts,gamma,proxfun);
    
    w = transpose(w);
    
    % normalize
    
    b = sqrt(sum(abs(w).^2,1))';
    inds_small = abs(b) < 10*eps(1)*max(b);
    b( inds_small ) = 1.0;
    w = w*diag(1./b);
    w(:,inds_small) = 0.0;
    
    if (nargout > 3)
        varargout{1} = w*diag(e)*pinv(w); % projected propagator
    end

    b( inds_small ) = 0.0;
    
    % unproject dmd modes
    
    w = u*w;
    
else
    
    % fit to all of data
    
    m = length(t);
    [is,~] = size(X);
    ia = r;
    n = r;
    
    [w,e,~,~,~,~] = varpro2(transpose(X),t, ...
        @varpro2expfun,@varpro2dexpfun,m,n,is,ia,alpha_init, ...
        opts,copts,gamma,proxfun);
    
    w = transpose(w);
    
    % normalize
    
    b = sqrt(sum(abs(w).^2,1))';
    inds_small = abs(b) < 10*eps(1)*max(b);
    b( inds_small ) = 1.0;    
    w = w*diag(1./b);
    w(:,inds_small) = 0.0;
    
    if (nargout > 3)
        wproj = u'*w;
        varargout{1} = wproj*diag(e)*pinv(wproj); %projected propagator
    end

    b( inds_small ) = 0.0;
    
    
end

if (nargout > 4)
    if (exist('u','var'))
        varargout{2} = u;
    else
        varargout{2} = [];
    end
end

if (nargout > 5)
    varargout{3} = w*diag(e)*pinv(w);
end



