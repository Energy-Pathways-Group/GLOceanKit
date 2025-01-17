function [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsWithFreeSurface(Finv,Ginv,h,Nj,Lz)
Nz = size(Finv,1);
nModes = size(Finv,2);

% Make these matrices invertible by adding the barotropic mode
% to F, and removing the lower boundary of G.
Finv = cat(2,ones(Nz,1),Finv); % [Nz Nj+1]
Ginv = Ginv(2:end,:);  % [Nz-1 Nj]

% Compute the precondition matrices (really, diagonals)
P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

% Now create the actual transformation matrices
PFinv = Finv./P;
QGinv = Ginv./Q;
PF = inv(PFinv); % [Nj+1 Nz]
QG = inv(QGinv); % [Nj Nz-1]

maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
if maxCond > 1000
    warning('Condition number is %f the vertical transformations.',maxCond);
end
% size(PFinv)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
% but, we will only multiply by vectors [Nj 1], so dump the
% last column. Now size(PFinv) = [Nz x Nj].
PFinv = PFinv(:,1:end-1);

% size(PF)=[Nj+1, Nz], but we don't care about the last mode
PF = PF(1:end-1,:);

% size(QGinv) = [Nz-1, Nj], need zeros for the lower boundaries
% and add the 0 barotropic mode, so size(G) = [Nz, Nj],
QGinv = QGinv(:,1:end-1); % dump Nyquist
QGinv = cat(2,zeros(Nz-1,1),QGinv); % add barotropic mode
QGinv = cat(1,zeros(1,Nj),QGinv); % add zeros at along the bottom

% Now have to do the same thing to the condition matrix
Q = cat(2,0,Q(1:end-1));

% size(QG) = [Nj, Nz-1], need a zero for the barotropic
% mode, but also need zeros for the boundary
QG = cat(1,zeros(1,self.Nz-1),QG(1:end-1,:)); % dump the Nyquist mode, add a barotropic mode (all zeros)
QG = cat(2,zeros(Nj,1), QG); % add the bottom boundary

% want size(h)=[1 1 Nj]
h = shiftdim(h,-1);

P = shiftdim(P(1:end-1),-1);
Q = shiftdim(Q,-1);
end