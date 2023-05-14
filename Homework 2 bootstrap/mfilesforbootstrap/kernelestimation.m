
function [yk additionalinfo] = kernelestimation(x,x0,varargin);

% kernelestimation - ThreshLab2/MultiscaleKernel/ --- Kernel smoothing routine
%   Usage
%      for kernel smoothing
%         yk = kernelestimation(x,x0,y0,h,kerneltype,side);
%         yk = kernelestimation(x,x0,y0,h,kerneltype,'reflex');
%      for kernel density estimation
%         yk = kernelestimation(x,X,h,kerneltype,side);   
%      for k-nearest neighbour estimation (Loftsgaarden & Quesenberry, 1965)
%         yk = kernelestimation(x,X,'knearest',k,'reflex');
%         yk = kernelestimation(x,x0,y0,'knearest',k,'reflex');
%      for adaptive kernel estimation (Breiman, Meisel, Purcell, 1977)
%         yk = kernelestimation(x,X,'adaptive',k,'reflex');
%         yk = kernelestimation(x,x0,y0,'adaptive',k,'reflex');
%   Inputs
%      x      x-coordinates in which kernel is evaluated
%      x0     x-coordinates of sample locations in smoothing problem
%      X      observations/sample vector in density estimation problem
%   Optional Inputs
%      y0     sample observations in smoothing problem
%      h      real scalar or vector; if h is a vector, then it must me
%             preceeded by vector y0, or, in the case of density estimation,
%             it must be preceeded by 'density'
%      'density'   sets density estimation flag to one (in most cases, the
%                  routine is able to detect whether the input is a density
%                  estimation problem without explicitly setting the flag to 1)
%      kerneltype  string-type variable, see evalkernel
%      side        -1,0,1,'L','l','Left','left','R','r','Right','right'
%      'reflex'   includes reflection near end points for bias correction
%      'k',k      number of neighbours in k-nearest neighbour estimation
%                 explicit declaration of 'k' is not always strictly necessary,
%                 but, for instance, if k is 1/2, then this value is identified
%                 with k only if k is preceeded by 'k'
%      'GM'       use Gasser-Müller estimation
%      'NW'       use Nadaraya-Watson estimation (default in smoothing)
%   Outputs
%      yk     kernel smoothing evaluated in x
%   Notes
%      reported bug: when the bandwidth is an integer, it may not be
%      interpreted well by the routine when its declaration is not preceeded by
%      'bandwidth'
%      x, x0 and y0 are supposed to be arranged as rows.
%      if y0 is a square matrix of size equal to the length of x0, then the
%      smoothing operation is performed on all ROWS.
%      As a consequence
%      kernelestimation(row(x),row(x),sparse(eye(N)),hj,kerneltype,side);
%      This version allows cut-off/one-sided kernels. (input variable side)
%   See also
%      help evalkernel

x = row(x);
x0 = row(x0);
additionalinfo = {};

nobs = length(x0);
neval = length(x);
% xx = [x(1) (x(1:neval-1)+x(2:neval))/2 x(neval)];
% dx = xx(2:neval+1)-xx(1:neval);


% [nrows ny] = size(y0);
% densityest = (nrows*ny<2);
% % if nrows == 1 & ny == 1, densityest = true; else, densityest = false; end
% if densityest == true, 
%    if nargin>4, varargin = {kerneltype,varargin{:}}; end
%    kerneltype = h;
%    h = y0;
%    y0 = ones(1,nobs)/nobs;
% else
%    if ny ~= nobs, 
%       if nrows ~= nobs, error('dimension mismatch x0 and y0')
%       else nrows = ny; y0 = y0'; end
%    end
% end
% if ischar(h),
%    varargin = {h,varargin{:}};
%    h = Inf;
% end

side = 0;
reflex = false;
nvarargin = length(varargin);
k = 1; adaptive = 0;
h = NaN;
arg = 1;
y0 = NaN;
hmin = 0;
tol = NaN;
densityest = false;
smoothing = 'NW';
makes0 = false;
while arg<=nvarargin, vararg = varargin{arg}; 
   if ischar(vararg), switch(vararg)
   case {'R','r','right','Right'},
      side = 1;
   case {'L','l','Left','left'},
      side = -1;
   case {'reflex','reflect'},
      reflex = true;
   case {'knearest'},
      k = max(k,1); adaptive = 1;
   case {'adaptive'},
      k = max(k,1); adaptive = 2;
   case {'knearest2','knearestLR'},
      k = max(k,1); adaptive = 3;
   case {'adaptive2','adaptiveLR'},
      k = max(k,1); adaptive = 4;
   case {'k'}, % in most cases, k is recognized without explicit declaration of
               % 'k'; explicit declaration is useful for noninteger values of k
      k = varargin{arg+1}; arg = arg+1;
   case {'bandwidth','h'},
      h = varargin{arg+1}; arg = arg+1;
   case {'hmin'},
      hmin = varargin{arg+1}; arg = arg+1;
   case {'tol'},
      tol = varargin{arg+1}; arg = arg+1;
      % when the kernel function is uniform, and k-nearest neighbours is used,
      % it may be usefull to add a small value (tol=eps) to the bandwidth in
      % the calculation of the neighbours covered by the kernel function
      % support (see below)
   case {'density','densityest','densityestimation'},
      densityest = true; smoothing = 'NW';
   case {'GM','Gasser-Muller','Gasser-Mueller','Gasser'},
      smoothing = 'GM';
      makes0 = true;
   case {'NW','Nadaraya-Watson'},
      smoothing = 'NW';
   otherwise
      kerneltype = vararg;
   end
   elseif isnumeric(vararg),
      if max(size(vararg)) == 1,
         r = vararg-floor(vararg);
         if abs(r)<eps 
            if adaptive > 0,
               k = max(k,vararg); 
            elseif abs(vararg)<2,
               side = vararg;
            elseif isnan(h(1,1)),
               h = vararg;
            end
         elseif isnan(h(1,1)),
            h = vararg;
         end
      else 
         if any(size(vararg)==nobs) & isnan(y0(1,1)) & densityest == false,
            y0 = vararg;
         elseif any(size(vararg)==nobs) | any(size(vararg)==neval),
            h = vararg;
         end
      end
   end
   arg = arg+1;
end
[nrows ny0] = size(y0);
if isnan(y0(1,1)), 
   densityest = true; 
else,
   [nrows ny0] = size(y0);
   if ny0 ~= nobs, 
      if nrows ~= nobs, error('dimension mismatch x0 and y0')
      else nrows = ny0; y0 = y0'; end
   end
end
if densityest == true, y0 = ones(1,nobs)/nobs; end
if makes0 == true,
   % a = min(x0); b = max(x0);
   a = -Inf; b = Inf;
   s0 = [a row(x0(2:end)+x0(1:end-1))/2 b];
end

if reflex==true,
   a = min(x); b = max(x);
   ll = find(x0-a<max(max(h))); rr = find(b-x0<max(max(h)));
   x0 = [2*a-x0(ll) x0 2*b-x0(rr)];
   y0 = [y0(ll) y0 y0(rr)];
   nobs = length(x0);
end

if isnan(tol),
   if ismember(kerneltype,{'uniform','Uniform','unif'}) && ...
      ismember(adaptive,[1 3]),
      tol = eps;
   else
      tol = 0;
   end
end

q = k-floor(k); k = floor(k); heval = true;
% The binary variable heval indicates whether the bandwidths in h are 
% associated to each given point in x or rather to each given point in x0. 
% In principle, the size of h reveals in which of the two situations we are, 
% unless x0 and x have the same size. This is why we provide this additional
% variable
switch(adaptive),
case 1,
   h = zeros(2,neval);
   for i=1:neval,
      di = sort(abs(x(i)-x0));
      h(1,i) = di(k)+q*(di(k+1)-di(k));
   end
   h(2,1:neval) = h(1,1:neval);
case 2,
   h = zeros(2,nobs); heval = false;
   for i=1:nobs,
      di = sort(abs(x0(i)-x0));
      h(1,i) = di(k+1)+q*(di(k+2)-di(k+1));
   end
   h(2,1:nobs) = h(1,1:nobs);
case 3,
   h = zeros(2,neval);
   for i=1:neval,
      di = x(i)-x0; di(di<0) = Inf; di = sort(di);
      h(1,i) = di(k)+q*(di(k+1)-di(k));
   end
   for i=1:neval,
      di = x0-x(i); di(di<0) = Inf; di = sort(di);
      h(2,i) = di(k)+(di(k+1)-di(k));
   end
case 4,
   h = zeros(2,nobs); heval = false;
   for i=2:nobs, % for i=1, the left bandwidth is not defined here; its value
                 % (zero for the moment) will be overwritten by
                 % h = max(h,hmin);
      di = x0(i)-x0; di(di<0) = Inf; di = sort(di);
      h(1,i) = di(k+1)+q*(di(k+2)-di(k+1));
   end
   for i=1:nobs-1, % for i=nobs, the right bandwidth is not defined here; its
                   % value (zero for the moment) will be overwritten by
                   % h = max(h,hmin);
      di = x0-x0(i); di(di<0) = Inf; di = sort(di);
      h(2,i) = di(k+1)+q*(di(k+2)-di(k+1));
   end
end

h = max(h,hmin);
hh = h;
if max(size(h)) == 1, hh = h*ones(2,neval); end
if max(size(h)) == 2,
   if min(size(h)) == 1, hh = [h(1)*ones(1,neval);h(2)*ones(1,neval)]; end
end
if ~ismember(size(hh,2),[nobs neval]), hh = hh'; end
if ~ismember(size(hh,2),[nobs neval]), 
   error('dimensions h and x/x0 mismatch'),
end
if size(hh,1) == 1, hh = [hh;hh]; end
nh = size(hh,2);

if side ==1, hh(1,1:nh) = 0; end
if side ==-1, hh(2,1:nh) = 0; end

kernelsupport = zeros(size(hh));
switch kerneltype
case {'Gaussian','gaussian','gauss','Gauss'}
   kernelsupport(abs(hh)>eps) = Inf;
otherwise
   kernelsupport = hh;
end

yk = zeros(nrows,neval);
sk = zeros(size(x));

for s = [-1,1],
   lr = (3+s)/2;
   for i = 1:nobs,
      % if nh == neval, kernelsupporti = kernelsupport(lr,1:nh);
      % elseif nh == nobs, kernelsupporti = kernelsupport(lr,i);
      % else, error('dimensions h and x/x0 mismatch'),
      if ismember(nh,[neval nobs])==false,
         error('dimensions h and x/x0 mismatch'),
      elseif heval == true, kernelsupporti = kernelsupport(lr,1:nh);
      else, kernelsupporti = kernelsupport(lr,i);
      end
      if makes0 == true,
         ii = (abs(x-s0(i))<kernelsupporti+tol | ...
               abs(x-s0(i+1))<kernelsupporti+tol);
      else
         ii = (abs(x-x0(i))<kernelsupporti+tol);
      end
      % we add tol because the bandwidth may be based on the distance to the
      % k nearest neighbours. In that case, we want all k nearest neighbours to
      % be taken into account, especially when the kernel function is the
      % uniform function.
      % keep points on one side only; be careful with points on the precise
      % spot of x0(i)
      % ii = ii & (s*(x-x0(i))>0);     this doesn't work for points very
      %                                close to x0(i)
      if s==1, ii = ii & ((x-x0(i))>-eps);
         else, ii = ii & ((x0(i)-x)>eps);
      end
      ii = find(ii);
      if ~isempty(ii),
         jj = find(abs(y0(1:nrows,i))>eps);
         if ismember(nh,[neval nobs])==false,
            error('dimensions h and x/x0 mismatch'),
         elseif heval == true,
            hiii = hh(lr,ii); h2iii = (hh(lr,ii)+hh(3-lr,ii))/2;
         else,
            hiii = hh(lr,i)*ones(size(ii));
            h2iii = (hh(lr,i)+hh(3-lr,i))*ones(size(ii))/2;
         end
         switch smoothing,
         case {'NW'}, % Nadaraya-Watson or density estimation
            Kx0ii = evalkernel(x(ii)-x0(i),hiii,kerneltype)./h2iii;
            yk(jj,ii) = yk(jj,ii)+y0(jj,i)*Kx0ii;
            sk(ii) = sk(ii)+Kx0ii;
         case {'GM'} % Gasser-Müller
            wii = integralkernel(x(ii)-s0(i),hiii,kerneltype)-...
                  integralkernel(x(ii)-s0(i+1),hiii,kerneltype);
            yk(jj,ii) = yk(jj,ii)+y0(jj,i)*wii;
            sk(ii) = 1;
         end
         % Ix0(i) = sum(Kx0ii.*dx(ii));
      end
   end
end
if densityest == false,
   % lines below: points with no prediction point within the bandwidth take the
   % closest observation as prediction value
   oo = find(abs(sk)<eps);
   if ~isempty(oo) & ...
      ~ismember(kerneltype,{'none','dirac','Dirac','kronecker','Kronecker'})
      disp(['points without prediction point within bandwidth:' num2str(oo)])
   end
   if ~isempty(oo), for j = row(oo),
      [xi i] = min(abs(x(j)-x0));
      yk(1:nrows,j) = y0(1:nrows,i);
      sk(j) = 1;
   end, end
   sk = repmat(sk,nrows,1);
   yk = yk./sk;
end
% additionalinfo = {hh};
