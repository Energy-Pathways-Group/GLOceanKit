% https://math.jhu.edu/~feilu/notes/DealiasingFFT.pdf
Lxyz = [2*pi 2*pi 1];
Nxyz = [16 16 2];
wvt = WVTransformConstantStratification(Lxyz, Nxyz,shouldAntialias=0);

[X,Y,Z] = wvt.xyzGrid;

Lx = wvt.Lx;

%%
antialiasMask = WVGeometryDoublyPeriodic.maskForAliasedModes(wvt.kAxis,wvt.lAxis,1);
trueAliasedModes = zeros(size(antialiasMask));
[K,L] = ndgrid(round(wvt.kAxis),round(wvt.lAxis));
validModes = find(antialiasMask==0);
for iIndex=1:length(validModes)
    k_i = K(validModes(iIndex));
    l_i = L(validModes(iIndex));
    for jIndex=iIndex:length(validModes)
        k_j = K(validModes(jIndex));
        l_j = L(validModes(jIndex));

        k_ij = k_i+k_j+floor(wvt.Nx/2);
        l_ij = l_i+l_j+floor(wvt.Ny/2);

        if k_ij < 0 || k_ij >= wvt.Nx ||  l_ij < 0 || l_ij >= wvt.Ny
            % this combination of modes aliases
            iAliasIndex = mod(k_ij,wvt.Nx) + 1;
            jAliasIndex = mod(l_ij,wvt.Ny) + 1;
            trueAliasedModes(iAliasIndex,jAliasIndex) = 1;
        end
    end
end

% Are there any modes that we flagged as being aliased, that were not
% flagged as being aliased?
any(antialiasMask(logical(trueAliasedModes))==0)

return

%% cos(alpha)*cos(beta) = 0.5*(cos(alpha+beta) + cos(alpha-beta))

% A correct, non-aliased implementation must have all combinations of
% quadratic multiplication
sin_resolved = zeros(1,wvt.Nkl);
sin_aliased = zeros(1,wvt.Nkl);
sin_zero = zeros(1,wvt.Nkl);
cos_resolved = zeros(1,wvt.Nkl);
cos_aliased = zeros(1,wvt.Nkl);
cos_zero = zeros(1,wvt.Nkl);

quadraticOperation{1} = @(sin_a,sin_b,cos_a,cos_b) sin_a.*cos_b + sin_b.*cos_a;
quadraticOperation{2} = @(sin_a,sin_b,cos_a,cos_b) sin_a.*cos_b - sin_b.*cos_a;
quadraticOperation{3} = @(sin_a,sin_b,cos_a,cos_b) cos_a.*cos_b - sin_a.*sin_b;
quadraticOperation{4} = @(sin_a,sin_b,cos_a,cos_b) cos_a.*cos_b + sin_a.*sin_b;

wavenumberOperation{1} = @(m_a,m_b) m_a + m_b;
wavenumberOperation{2} = @(m_a,m_b) m_a - m_b;
wavenumberOperation{3} = @(m_a,m_b) m_a + m_b;
wavenumberOperation{4} = @(m_a,m_b) m_a - m_b;

op_resolved = zeros(1,wvt.Nkl);
op_aliased = zeros(1,wvt.Nkl);
op_zero = zeros(1,wvt.Nkl);

for m=1:wvt.Nkl
    for n=1:wvt.Nkl
        k_m = wvt.k(m); l_m = wvt.l(m);
        k_n = wvt.k(n); l_n = wvt.l(n);
        sin_m = sin(k_m*X + l_m*Y);
        cos_m = cos(k_m*X + l_m*Y);
        sin_n = sin(k_n*X + l_n*Y);
        cos_n = cos(k_n*X + l_n*Y);

        for iOp=1:4
            quad_mn = quadraticOperation{iOp}(sin_m,sin_n,cos_m,cos_n);
            quad_mn_bar = wvt.transformFromSpatialDomainWithFourier(quad_mn);
            quad_nm_index = find(abs(quad_mn_bar(1,:)) > 1e-6);
            if isempty(quad_nm_index)
                op_zero(m) = op_zero(m) + 1;
            elseif length(quad_nm_index) > 1
                error('Returned two wavenumbers!')
            else
                k_mn = wvt.k(quad_nm_index); l_mn = wvt.l(quad_nm_index);
                if k_mn == wavenumberOperation{iOp}(k_m,k_n) && l_mn == wavenumberOperation{iOp}(l_m,l_n)
                    op_resolved(m) = op_resolved(m) + 1;
                else
                    op_aliased(m) = op_aliased(m) + 1;
                end
            end
        end

        % sin_mn = sin_m.*cos_n+sin_n.*cos_m;
        % cos_mn = cos_m.*cos_n-sin_n.*sin_m;
        % 
        % sin_mn_bar = wvt.transformFromSpatialDomainWithFourier(sin_mn);
        % cos_mn_bar = wvt.transformFromSpatialDomainWithFourier(cos_mn);
        % 
        % sin_nm_index = find(abs(sin_mn_bar(1,:)) > 1e-6);
        % cos_nm_index = find(abs(cos_mn_bar(1,:)) > 1e-6);
        % 
        % if isempty(sin_nm_index)
        %     sin_zero(m) = sin_zero(m) + 1;
        % elseif length(sin_nm_index) > 1
        %     error('Returned two wavenumbers!')
        % else
        %     k_mn = wvt.k(sin_nm_index); l_mn = wvt.l(sin_nm_index);
        %     if k_mn == (k_m + k_n) && l_mn == (l_m + l_n)
        %         sin_resolved(m) = sin_resolved(m) + 1;
        %     else
        %         sin_aliased(m) = sin_aliased(m) + 1;
        %     end
        % end
        % 
        % if isempty(cos_nm_index)
        %     cos_zero(m) = cos_zero(m) + 1;
        % elseif length(cos_nm_index) > 1
        %     error('Returned two wavenumbers!')
        % else
        %     k_mn = wvt.k(cos_nm_index); l_mn = wvt.l(cos_nm_index);
        %     if k_mn == (k_m + k_n) && l_mn == (l_m + l_n)
        %         cos_resolved(m) = cos_resolved(m) + 1;
        %     else
        %         cos_aliased(m) = cos_aliased(m) + 1;
        %     end
        % end
    end
end

%%
op_total = op_aliased + op_resolved + op_zero;
[op_resolved_kl,op_zero_kl,op_aliased_kl,op_total_kl] = wvt.transformToKLAxes(op_resolved,op_zero,op_aliased,op_total);
figure
tiledlayout(1,3)
nexttile
jpcolor(wvt.kAxis,wvt.lAxis,op_resolved_kl.'), shading flat
nexttile
jpcolor(wvt.kAxis,wvt.lAxis,op_zero_kl.'), shading flat
nexttile
jpcolor(wvt.kAxis,wvt.lAxis,op_aliased_kl.'), shading flat

%%
k_n = 3;
kx = 2*pi*k_n/Lx;
f = cos(kx*X);
f2 = 0.5 + 0.5*cos(2*kx*X);
f2_actual = wvt.transformFromSpatialDomainWithFourier(f.*f);
f2_actual_2d = wvt.transformToKLAxes(f2_actual);

% figure, plot(f2(1,:))

%%
n=5; m=5; u2_bar = fft(fft(cos(n*wvt.dk*X).*cos(m*wvt.dk*X),wvt.Nx,1),wvt.Ny,2)/wvt.Nx/wvt.Ny;
n=5; m=5; u2_bar = fft(fft(0.5*(cos((n+m)*wvt.dk*X) + cos((n-m)*wvt.dk*X)),wvt.Nx,1),wvt.Ny,2)/wvt.Nx/wvt.Ny;
p = mod(n+m,wvt.Nx);
p_index = p+1;
