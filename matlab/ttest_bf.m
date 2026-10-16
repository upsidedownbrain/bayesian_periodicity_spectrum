function [bf10, pValue, T, N, log10bf, diagnostics] = ttest_bf(X, varargin)
% TTEST_BF Bayesian t-test using a noncentral-t likelihood.
%
% Authors: Jesús Pardo-Valencia, Guglielmo Foffani
% Version: 2.0.0 
% Date: October 2026
% Repository: GitHub
% License: MIT License
% Copyright (c) 2025 jpv
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
% Description:
%   Computes BF10 for one-sample, paired, or independent-samples t-tests.
%
%   Under H1:
%
%       T | delta ~ noncentral-t(df, sqrt(N)*delta)
%
%   Priors:
%       tail = 'right' : normalized half-Cauchy, delta > 0
%       tail = 'left'  : normalized half-Cauchy, delta < 0
%       tail = 'both'  : full Cauchy
%
%   Default Cauchy scale:
%
%       r = sqrt(2)/2
%
%   The exact JZS normal/inverse-gamma mixture is integrated in log(g).
%   Cached composite quadrature is checked at two resolutions and uses
%   analytic bounds for omitted tails. Difficult cases use adaptive
%   quadrature. Failure to meet RelTol raises PrecisionNotAchieved.
%   This formulation avoids nctpdf without changing the Cauchy prior.
%
% Optional ESS correction:
%   For one-sample or paired analyses:
%
%       Neff = N*(1-rho1)/(1+rho1)
%
%       T_eff = T*sqrt(Neff/N)
%
%   This is an AR(1)-based approximation.
%
% Usage:
%
%   One-sample:
%       [bf10,p,T,N] = ttest_bf(X)
%
%   One-sample against M:
%       [bf10,p,T,N] = ttest_bf(X,M)
%
%   Paired:
%       [bf10,p,T,N] = ttest_bf(X,Y,'paired',true)
%
%   Independent samples:
%       [bf10,p,T,N] = ttest_bf(X,Y,'paired',false)
%
%   Direct T (one-sample/paired; rho1 applies ESS to nominal T and N):
%       [bf10,p,T,N] = ttest_bf('T',Tvalue,'N',Nvalue)
%
% Optional parameters:
%   'tail'   : 'both', 'right', 'left'
%   'scale'  : Cauchy scale
%   'T'      : supplied t statistic(s)
%   'N'      : sample size or effective sample size
%   'rho1'   : lag-1 autocorrelation, scalar or one value per test
%   'Dim'    : sample dimension (1 or 2); default first nonsingleton
%              dimension, with vector samples normalized to columns
%   'paired' : true/false
%   'nQuad'  : initial nodes per log-g panel (default 16; checked/refined)
%              Its meaning differs from the old theta-domain quadrature.
%   'RelTol' : target relative integration accuracy (default 1e-10)
%   'MaxQuad': maximum nodes per panel before adaptive fallback
%              (default max(128,2*nQuad))
%
%   Direct independent T: use scalar T and 'N',[n1 n2].
%   Independent tests use pooled variance and do not accept rho1.
%   NaN observations are omitted; paired tests omit incomplete pairs.
%   Non-finite t statistics are rejected rather than returning invalid BF.
%
% Outputs:
%   bf10    : Bayes factor H1/H0
%   pValue  : frequentist p-value
%   T       : t statistic used
%   N       : valid N or Neff per test; independent tests return n1*n2/(n1+n2)
%
%   log10bf : log10(BF), retained even when raw BF over/underflows
%   diagnostics : converged, relativeErrorEstimate, nodesPerPanel,
%                 adaptiveFallback, bfOverflow, bfUnderflow
%   The error estimate checks quadrature and omitted tails; it is not
%   a validation of the ESS approximation or other statistical assumptions.
%


%% ========================================================================
% Parse positional inputs
% ========================================================================

if isnumeric(X)

    if mod(numel(varargin),2) == 0
        Y = 0;
        parms = varargin;
    else
        Y = varargin{1};

        if numel(varargin) > 1
            parms = varargin(2:end);
        else
            parms = {};
        end
    end

else

    % Example:
    % ttest_bf('T',[1 2 3],'N',30,'tail','right')

    parms = [{X}, varargin];
    X = [];
    Y = [];

end


%% ========================================================================
% Parse parameters
% ========================================================================

p = inputParser;

p.addParameter( ...
    'tail','both', ...
    @(x) (ischar(x) || isstring(x)) && ...
    ismember(lower(string(x)),["both","right","left"]) );

p.addParameter( ...
    'scale',sqrt(2)/2, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x > 0 );

p.addParameter( ...
    'T',[], ...
    @isnumeric );

p.addParameter( ...
    'N',[], ...
    @isnumeric );

p.addParameter( ...
    'rho1',[], ...
    @(x) isempty(x) || (isnumeric(x) && isreal(x) && ...
    all(isfinite(x(:))) && all(abs(x(:)) < 1)) );

p.addParameter( ...
    'paired',true, ...
    @(x) islogical(x) && isscalar(x) );

p.addParameter( ...
    'nQuad',16, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 8 && x <= 128 && x == fix(x) );

p.addParameter('RelTol',1e-10, ...
    @(x) isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 1e-12 && x <= 1e-3);
p.addParameter('MaxQuad',[], ...
    @(x) isempty(x) || (isnumeric(x) && isreal(x) && isscalar(x) && isfinite(x) && x >= 16 && x <= 256 && x == fix(x)));

p.addParameter('Dim',[], ...
    @(x) isempty(x) || (isnumeric(x) && isscalar(x) && ismember(x,[1 2])));

p.parse(parms{:});

tail   = char(lower(string(p.Results.tail)));
r      = p.Results.scale;
rho1   = p.Results.rho1;
paired = p.Results.paired;
nQuad = p.Results.nQuad;
relTol = p.Results.RelTol;
maxQuad = p.Results.MaxQuad;
if isempty(maxQuad),maxQuad=max(128,2*nQuad);end
if maxQuad<=nQuad,error('ttest_bf:InvalidMaxQuad','MaxQuad must exceed nQuad to verify convergence.');end


%% ========================================================================
% Obtain T, N, df
% ========================================================================

if isempty(p.Results.T)
    % Vectors are samples; matrices are separate tests along Dim.
    validateattributes(X,{'numeric'},{'real','nonempty','2d'},mfilename,'X');
    validateattributes(Y,{'numeric'},{'real','nonempty','2d'},mfilename,'Y');
    if any(isinf(X(:))) || any(isinf(Y(:)))
        error('ttest_bf:InfiniteData','Data must not contain Inf. NaN denotes missing observations.');
    end
    if ~isempty(p.Results.N)
        error('ttest_bf:UnusedN','N is calculated from the data; supply N only with direct T input.');
    end
    if isvector(X) && (isscalar(Y) || isvector(Y)) && isempty(p.Results.Dim)
        X = X(:);
        if ~isscalar(Y), Y = Y(:); end
    end
    dim = p.Results.Dim;
    if isempty(dim)
        dim = find(size(X) ~= 1,1);
        if isempty(dim), dim = 1; end
    end

    if isscalar(Y)
        if ~isfinite(Y)
            error('ttest_bf:InvalidNull','The one-sample null mean must be finite.');
        end
        N = sum(~isnan(X),dim);
        if any(N(:) < 2)
            error('ttest_bf:InsufficientData','Each test requires at least two valid observations.');
        end
        [~,~,~,stats] = ttest(X,Y,'Dim',dim);
        useESS = true;
    elseif paired
        if ~isequal(size(X),size(Y))
            error('ttest_bf:PairedSize','Paired matrices must have identical sizes and aligned observations.');
        end
        N = sum(~isnan(X) & ~isnan(Y),dim);
        if any(N(:) < 2)
            error('ttest_bf:InsufficientData','Each paired test requires at least two valid pairs.');
        end
        [~,~,~,stats] = ttest(X,Y,'Dim',dim);
        useESS = true;
    else
        otherDim = 3-dim;
        if size(X,otherDim) ~= size(Y,otherDim)
            error('ttest_bf:IndependentSize','Independent matrices must match outside the sample dimension.');
        end
        n1 = sum(~isnan(X),dim);
        n2 = sum(~isnan(Y),dim);
        if any(n1(:) < 2) || any(n2(:) < 2)
            error('ttest_bf:InsufficientData','Each independent sample requires at least two valid observations.');
        end
        [~,~,~,stats] = ttest2(X,Y,'Dim',dim);
        N = (n1 .* n2) ./ (n1 + n2);
        useESS = false;
    end
    T = stats.tstat;
    df = stats.df;
else
    % Direct T uses nominal N when rho1 is supplied; otherwise N may
    % already be effective. Do not supply rho1 again after correcting T/N.
    if ~isempty(X)
        error('ttest_bf:MixedInput','Use data input or direct T input, not both.');
    end
    if ~isempty(p.Results.Dim)
        error('ttest_bf:UnusedDim','Dim applies only to data input.');
    end
    T = p.Results.T;
    N = p.Results.N;
    validateattributes(T,{'numeric'},{'real','finite','nonempty','2d'},mfilename,'T');
    validateattributes(N,{'numeric'},{'real','finite','nonempty','positive'},mfilename,'N');
    if numel(N) == 2 && isscalar(T)
        % Legacy direct independent-samples form: N = [n1 n2].
        if any(N(:) < 2) || any(N(:) ~= fix(N(:)))
            error('ttest_bf:InvalidSampleSizes','Independent sample sizes must be integers >= 2.');
        end
        n1 = N(1); n2 = N(2);
        df = n1+n2-2;
        N = n1*n2/(n1+n2);
        useESS = false;
    else
        if any(N(:) <= 1)
            error('ttest_bf:InvalidN','One-sample or paired N must be greater than one.');
        end
        df = N-1;
        useESS = true;
    end
end

% Match parameter shapes before tcdf or element-wise ESS calculations.
originalSize = size(T);
if any(~isfinite(T(:)))
    error('ttest_bf:DegenerateData','The t statistic is not finite; check sample variance and valid observations.');
end
if ~isscalar(N)
    if numel(N) ~= numel(T)
        error('ttest_bf:ParameterSize','N must be scalar or have the same number of elements as T.');
    end
    N = reshape(N,originalSize);
end
if ~isscalar(df), df = reshape(df,originalSize); end

if ~isempty(rho1)
    if ~useESS
        error('ttest_bf:UnsupportedESS','rho1 is supported only for one-sample or paired tests.');
    end
    if ~isscalar(rho1)
        if numel(rho1) ~= numel(T)
            error('ttest_bf:ParameterSize','rho1 must be scalar or have the same number of elements as T.');
        end
        rho1 = reshape(rho1,originalSize);
    end
    if any(N(:)<2)
        error('ttest_bf:InvalidNominalN','ESS requires nominal N >= 2.');
    end
    r1 = max(min(rho1,0.99),-0.99);
    N_original = N;
    Neff = max(2,min(N_original,N_original .* (1-r1) ./ (1+r1)));
    T = T .* sqrt(Neff ./ N_original);
    N = Neff;
    df = Neff-1;
end

switch tail
    case 'both'
        pValue = 2 .* (1-tcdf(abs(T),df));
    case 'right'
        pValue = 1-tcdf(T,df);
    case 'left'
        pValue = tcdf(T,df);
end

T = T(:);
nT = numel(T);
if isscalar(N), N_vec = repmat(N,nT,1); else, N_vec = N(:); end
if isscalar(df), df_vec = repmat(df,nT,1); else, df_vec = df(:); end


% Log-domain JZS integration with convergence and truncation checks.
[logBF, diagnostics] = jzs_backend(T,N_vec,df_vec,r,tail,nQuad,relTol,maxQuad);
bf10 = exp(logBF);
log10bf = reshape(logBF/log(10),originalSize);
fields = fieldnames(diagnostics);
for k=1:numel(fields)
    diagnostics.(fields{k})=reshape(diagnostics.(fields{k}),originalSize);
end

%% ========================================================================
% Restore original shape
% ========================================================================

bf10   = reshape(bf10,originalSize);
pValue = reshape(pValue,originalSize);
T      = reshape(T,originalSize);

end



function [logBF,diag] = jzs_backend(t,N,nu,r,tail,q0,tol,qmax)
% Exact Cauchy/normal-IG representation of the JZS Bayes factor.
% Composite quadrature is checked at two resolutions, with explicit
% bounds on both omitted tails. Unresolved cases use adaptive integration.
t=t(:);N=N(:);nu=nu(:);n=numel(t);
logBF=nan(n,1);err=inf(n,1);used=zeros(n,1);fallback=false(n,1);
L=-12;U=max(48,2*ceil(max(2*log(max(abs(t),realmin))-log(N)-2*log(r)+40)/2));
if U-L>2048
    pending=(1:n)';
else
    pending=(1:n)';
    for rangeAttempt=1:8
        prev=quad_log(t(pending),N(pending),nu(pending),r,tail,q0,L,U);
        q=q0;
        while q<qmax && ~isempty(pending)
            q=min(2*q,qmax);
            cur=quad_log(t(pending),N(pending),nu(pending),r,tail,q,L,U);
            e=abs(expm1(prev-cur));
            [lb,ub]=tail_bounds(t(pending),N(pending),nu(pending),r,tail,L,U);
            trunc=exp(lb-cur)+exp(ub-cur);
            ok=isfinite(cur) & e<=tol/4 & trunc<=tol/4;
            logBF(pending(ok))=cur(ok);err(pending(ok))=e(ok)+trunc(ok);used(pending(ok))=q;
            pending=pending(~ok);prev=cur(~ok);
        end
        if isempty(pending),break;end
        % Expand only for a truncation problem; lack of resolution goes
        % directly to the independently refined adaptive calculation.
        [lb,ub]=tail_bounds(t(pending),N(pending),nu(pending),r,tail,L,U);
        if all(exp(lb-prev)+exp(ub-prev)<=tol/4),break;end
        L=L-2;U=U+16;
        if U-L>2048,break;end
    end
end
for j=pending.'
    [logBF(j),err(j)]=adaptive_log(t(j),N(j),nu(j),r,tail,tol);
    fallback(j)=true;
end
if any(~isfinite(logBF)) || any(err>tol)
    error('ttest_bf:PrecisionNotAchieved','BF integration did not meet RelTol; no unverified value is returned.');
end
diag=struct('converged',true(n,1),'relativeErrorEstimate',err, ...
    'nodesPerPanel',used,'adaptiveFallback',fallback, ...
    'bfOverflow',logBF>log(realmax),'bfUnderflow',logBF<log(realmin));
end

function y=quad_log(t,N,nu,r,tail,q,L,U)
[x,w]=legendre(q);edges=L:2:U;
v=reshape(x+edges(1:end-1)+1,1,[]);
weights=repmat(w,numel(edges)-1,1).';
y=zeros(numel(t),1);
chunk=max(1,floor(1e6/numel(v)));
for k=1:chunk:numel(t)
    ix=k:min(k+chunk-1,numel(t));
    z=log_integrand(v,t(ix),N(ix),nu(ix),r,tail)+log(weights);
    peak=max(z,[],2);
    y(ix)=peak+log(sum(exp(z-peak),2));
end
end

function y=log_integrand(v,t,N,nu,r,tail)
loga=log(N)+2*log(r)+v;
log1pa=softplus(loga);
logratio=2*log(abs(t))-log(nu);
logD0=softplus(logratio);
logDg=softplus(logratio-log1pa);
logx=logDg-logD0;
alpha=(nu+1)/2;
base=-.5*log1pa-alpha.*logx;
if ~strcmp(tail,'both')
    if strcmp(tail,'left'),t=-t;end
    logz=log(abs(t))+.5*(loga-log1pa)+.5*(log1p(1./nu)-logDg);
    z=exp(min(logz,log(realmax)));
    if all(nu==nu(1)),cdfdf=nu(1)+1;
    else,cdfdf=repmat(nu+1,1,numel(v));end
    adverse=t<0;
    signedZ=z.*sign(t);
    prob=tcdf(signedZ,cdfdf);
    ydir=base+log(2)+log(prob);
    % In the adverse tail, combine the incomplete beta with x^(-alpha)
    % analytically. This cancels the large powers before rounding, and
    % avoids underflow of tcdf for large negative z.
    mask=(adverse+zeros(size(v))) & (-alpha.*logx>32 | prob==0);
    if any(mask(:))
        aa=alpha+zeros(size(v));
        xx=exp(logx(mask));
        if any(xx==0)
            % x->0: I_x(a,1/2)/x^a -> 1/(a*B(a,1/2)).
            cf=zeros(size(xx));nz=xx>0;
            av=aa(mask);cf(nz)=log_beta_cf(av(nz),xx(nz));
        else
            cf=log_beta_cf(aa(mask),xx);
        end
        ydir(mask)=-.5*log1pa(mask)+.5*log1p(-xx) ...
            -log_beta_half(aa(mask))-log(aa(mask))+cf;
    end
    base=ydir;
end
y=base-.5*log(2*pi)-.5*v-.5*exp(-v);
end

function [lo,hi]=tail_bounds(t,N,nu,r,tail,L,U)
logD0=softplus(2*log(abs(t))-log(nu));
C=(nu+1)/2.*logD0;
factor=0;
if ~strcmp(tail,'both')
    factor=log(2);
    if strcmp(tail,'left'),t=-t;end
    C(t<=0)=0;
end
% For an adverse directional test BF(g)<=1/sqrt(1+N*r^2*g).
% For the other tests replace it by (1+t²/nu)^((nu+1)/2)
% times that bound. Prior density in log(g) <= exp(-v/2)/sqrt(2*pi).
hi=C-.5*(log(N)+2*log(r))-.5*log(2*pi)-U+factor;
% Mills bound for P(g<exp(L)), g~InvGamma(1/2,1/2).
lo=C+.5*L-.5*exp(-L)+.5*log(2/pi)+factor;
end

function [value,errorEstimate]=adaptive_log(t,N,nu,r,tail,tol)
L=-12;U=max(48,2*ceil((2*log(max(abs(t),realmin))-log(N)-2*log(r)+40)/2));
fun=@(v) log_integrand(v,t,N,nu,r,tail);
for k=1:12
    % A broad scan brackets the mode even when the scale/effect moves it.
    scan=linspace(L,U,max(61,ceil((U-L)/2)+1));ys=fun(scan);[~,i]=max(ys);
    a=scan(max(1,i-1));b=scan(min(numel(scan),i+1));
    [mode,np]=fminbnd(@(v)-fun(v),a,b,optimset('TolX',1e-9));peak=-np;
    wp=unique([mode mode+(-8:2:8)]);wp=wp(wp>L & wp<U);
    f=@(v)exp(fun(v)-peak);
    [v1,e1]=quadgk(f,L,U,'Waypoints',wp,'RelTol',max(tol/32,100*eps),'AbsTol',tol/32);
    [v2,e2]=quadgk(f,L,U,'Waypoints',unique([wp (L+U)/2]),'RelTol',max(tol/128,100*eps),'AbsTol',tol/128);
    value=peak+log(v2);
    [lb,ub]=tail_bounds(t,N,nu,r,tail,L,U);
    trunc=exp(lb-value)+exp(ub-value);
    errorEstimate=abs(v1/v2-1)+max(e1,e2)/v2+trunc;
    if isfinite(value) && errorEstimate<=tol,return;end
    L=L-2;U=U+16;
end
error('ttest_bf:PrecisionNotAchieved','Adaptive JZS integration did not converge.');
end

function y=softplus(x)
y=max(x,0)+log1p(exp(-abs(x)));
end

function [x,w]=legendre(n)
persistent counts xs ws
if isempty(counts),counts=[];xs={};ws={};end
i=find(counts==n,1);
if ~isempty(i),x=xs{i};w=ws{i};return;end
k=(1:n-1)';beta=k./sqrt(4*k.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));[x,ix]=sort(diag(D));w=2*(V(1,ix).^2)';
counts(end+1)=n;xs{end+1}=x;ws{end+1}=w;
end

function logh=log_beta_cf(a,x)
% Continued fraction for I_x(a,1/2), evaluated without its small prefactor.
b=.5;qab=a+b;qap=a+1;qam=a-1;
c=ones(size(a));d=1-qab.*x./qap;d=clamp(d);d=1./d;h=d;
active=true(size(a));
for m=1:4096
    ix=find(active);if isempty(ix),break;end
    ai=a(ix);xi=x(ix);m2=2*m;
    aa=m*(b-m).*xi./((qam(ix)+m2).*(ai+m2));
    di=clamp(1+aa.*d(ix));ci=clamp(1+aa./c(ix));di=1./di;hi=h(ix).*di.*ci;
    aa=-(ai+m).*(qab(ix)+m).*xi./((ai+m2).*(qap(ix)+m2));
    di=clamp(1+aa.*di);ci=clamp(1+aa./ci);di=1./di;delta=di.*ci;
    h(ix)=hi.*delta;d(ix)=di;c(ix)=ci;
    active(ix(abs(delta-1)<5e-14))=false;
end
if any(active) || any(h<=0) || any(~isfinite(h))
    error('ttest_bf:PrecisionNotAchieved','Directional incomplete-beta calculation did not converge.');
end
logh=log(h);
end
function y=clamp(y)
small=abs(y)<1e-300;y(small)=1e-300.*(1-2*(y(small)<0));
end

function y=log_beta_half(a)
% Avoid subtracting large gammaln values for B(a,1/2).
y=betaln(a,.5);large=a>100;x=a(large);
y(large)=.5*log(pi)-.5*log(x)+1./(8*x)-1./(192*x.^3)+1./(640*x.^5);
end
