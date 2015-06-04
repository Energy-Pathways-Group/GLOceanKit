function [varargout] = SecondMomentMatrix( varargin )
%% SecondMomentMatrix
%
%	Jeffrey J. Early, 2014
%
% Compute the components of the second moment for a group of particles.
%
% Calling
%	[M_qq, M_rr, M_qr] = SecondMomentMatrix( x, y )
% will return the matrix components in center of mass coordinates, (q,r).
%
% Calling
%	[minD, maxD, theta] = SecondMomentMatrix( x, y, 'eigen' )
% will return the matrix components in center of mass coordinates, but in the eigenbasis.
%
% Note that if you want to compute diffusivity by assuming isotrop, minD+maxD gives you
% the dispersion along each dimension, so you'll want 1/4 the slope to get kappa.

x=varargin{1};
y=varargin{2};
if nargin > 2
    basis = varargin{3};
else
    basis = '';
end

[~, ~, q, r] = CenterOfMass( x, y );

% Create the second moment matrix
M_qq = mean(q.*q,2);
M_rr = mean(r.*r,2);
M_qr = mean(q.*r,2);

if ( strcmp(basis, 'eigen') || strcmp(basis, 'eigenbasis') )
    nT = size(x,1);
	theta = zeros(nT,1);
	minD = zeros(nT,1);
	maxD = zeros(nT,1);
	for iTime =1:nT
		[A, B] = eig([M_qq(iTime), M_qr(iTime); M_qr(iTime), M_rr(iTime)]);
		v = A(:,2);
		theta(iTime) = atan(v(2)/v(1));
		minD(iTime) = B(1,1);
		maxD(iTime) = B(end,end);
    end
	
    varargout{1} = minD;
    varargout{2} = maxD;
    varargout{3} = theta;
	return
else
    varargout{1} = M_qq;
    varargout{2} = M_rr;
    varargout{3} = M_qr;
	return
end