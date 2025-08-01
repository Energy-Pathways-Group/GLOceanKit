function [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,Lz)
% return the normalized projection operators with prefactors
%
% This function uses InternalModesWKBSpectral to compute the
% quadrature points of a given stratification profile.
%
% - Topic: Internal
% - Declaration: [P,Q,PFinv,PF,QGinv,QG,h,w] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,Lz)
% - Returns P: preconditioner for F, size(P)=[Nj 1]
% - Returns Q: preconditioner for Q, size(Q)=[Nj 1]
% - Returns PFinv: normalized Finv, size(PFinv)=[Nz x Nj]
% - Returns PF: normalized F, size(PF)=[Nj x Nz]
% - Returns QGinv: normalized QGinv, size(QGinv)=[Nz x Nj]
% - Returns QG: normalized G, size(QG)=[Nj x Nz]
% - Returns h: eigenvalue, size(h)=[Nj x 1]
% - Returns w: eigenvalue, size(h)=[Nj x 1]
Nz = size(Finv,1);
nModes = size(Finv,2);

% Make these matrices invertible by adding the barotropic mode
% to F, and removing the boundaries of G.
Finv = cat(2,ones(Nz,1),Finv);
Ginv = Ginv(2:end-1,1:end-1);

% Compute the precondition matrices (really, diagonals)
P = max(abs(Finv),[],1); % ones(1,size(Finv,1)); %
Q = max(abs(Ginv),[],1); % ones(1,size(Ginv,1)); %

% Now create the actual transformation matrices
PFinv = Finv./P;
QGinv = Ginv./Q;
PF = inv(PFinv);
QG = inv(QGinv);

b = zeros(Nz,1);
b(1) = Lz;
w = (PFinv.')\b;

maxCond = max([cond(PFinv), cond(QGinv), cond(PF), cond(QG)],[],2);
if maxCond > 1000
    warning('Condition number is %f the vertical transformations.',maxCond);
end
% size(F)=[Nz x Nj+1], barotropic mode AND extra Nyquist mode
% but, we will only multiply by vectors [Nj 1], so dump the
% last column. Now size(Fp) = [Nz x Nj].
PFinv = PFinv(:,1:end-1);

% size(Finv)=[Nj+1, Nz], but we don't care about the last mode
PF = PF(1:end-1,:);

% size(G) = [Nz-2, Nj-1], need zeros for the boundaries
% and add the 0 barotropic mode, so size(G) = [Nz, Nj],
QGinv = cat(2,zeros(Nz,1),cat(1,zeros(1,nModes-1),QGinv,zeros(1,nModes-1)));

% size(Ginv) = [Nj-1, Nz-2], need a zero for the barotropic
% mode, but also need zeros for the boundary
QG = cat(2,zeros(nModes,1), cat(1,zeros(1,Nz-2),QG),zeros(nModes,1));

% want size(h)=[Nj 1]
h = cat(1,1,reshape(h(1:end-1),[],1)); % remove the extra mode at the end

P = reshape(P(1:end-1),[],1);
Q = reshape(cat(2,1,Q),[],1);

PFinv = PFinv(:,1:Nj);
PF = PF(1:Nj,:);
P = P(1:Nj,1);
QGinv = QGinv(:,1:Nj);
QG = QG(1:Nj,:);
Q = Q(1:Nj,1);
h = h(1:Nj,1);
end