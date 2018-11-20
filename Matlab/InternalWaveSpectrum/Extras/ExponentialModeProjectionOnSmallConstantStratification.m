% Project the exponential modes onto the constant stratification modes for
% some chunk of the water column. It turns out, if you choose the upper 650
% meters (b/2), then, by several metrics, you produce a situation where the
% n_exp-th mode projects onto the n_small/2 constant mode.
%
% So, I think that if you then use a domain with 1300m depth then because
% the n-th mode in 1300 m projects onto the n_big*pi/1300 = n_small*pi/650,
% which means that 2*n_small = n_big so now n_exp = n_big, just like wkb
% scaling might suggest.

Lz = 4000;
z = linspace(-Lz,0,1000)';
im = InternalModesExponentialStratification([5.2e-3 1300],[-Lz 0],z,33);
[F,G,h] = im.ModesAtWavenumber(0);

z_const = linspace(-1300,0,1000)';
imConst = InternalModesConstantStratification(5.2e-3,[-1300 0],z_const,33);
[Fconst,Gconst,hconst] = imConst.ModesAtWavenumber(0);

return

upperIndices = find(z>-650);

modes = (1:20)';

weighted = zeros(length(modes),1);
for j=1:length(modes)
    [xbar,f] = CosineTransformForward(z(upperIndices),F(upperIndices,modes(j)));
    weighted(j) = sum(abs(xbar(modes)).*(modes))/sum(abs(xbar(modes)));
    weighted(j) = sum((xbar(modes).^2).*(modes))/sum(abs((xbar(modes).^2)));
    
%     [m,weighted(j)] =  max(abs(xbar));
end

figure, plot(weighted)
xlabel('modes')
ylabel('weighted projection')

[xbar,f] = CosineTransformForward(z(upperIndices),F(upperIndices,1));
figure
plot(abs(xbar(1:14))), ylog, hold on
for j=2:6
    [xbar,f] = CosineTransformForward(z(upperIndices),F(upperIndices,j));
    plot(abs(xbar(1:14)))
end
legend('1','2','3','4','5','6')