file = '/Users/jearly/Documents/Data/ForcedDissipativeQGTurbulence/lowResSSH.nc';
k = ncread(file, 'k');
l = ncread(file, 'l');
f = ncread(file, 'ssh_realp') + sqrt(-1)*ncread(file, 'ssh_imagp');
ssh1 = f;
l = cat(1,l(1:end-1), -flipdim(l(2:end,:),1));
ssh1Matlab = cat(1,f(1:end-1,:), conj(cat(2, flipdim(f(2:end,1),1), flipdim(flipdim(f(2:end,2:end),1),2))) )*length(k)*length(l);
ssh1t = ifft2(ssh1Matlab,'symmetric');

sum(sum(ssh1Matlab.*conj(ssh1Matlab)))/(length(k)*length(l))
sum(sum(ssh1t.*ssh1t))

file = '/Users/jearly/Documents/Data/ForcedDissipativeQGTurbulence/highResSSH.nc';
k = ncread(file, 'k');
l = ncread(file, 'l');
f = ncread(file, 'ssh_realp') + sqrt(-1)*ncread(file, 'ssh_imagp');
ssh2 = f;
l = cat(1,l(1:end-1), -flipdim(l(2:end,:),1));
ssh2Matlab = cat(1,f(1:end-1,:), conj(cat(2, flipdim(f(2:end,1),1), flipdim(flipdim(f(2:end,2:end),1),2))) )*length(k)*length(l);
ssh2t = ifft2(ssh2Matlab,'symmetric');

sum(sum(ssh2Matlab.*conj(ssh2Matlab)))/(4*length(k)*length(l))
sum(sum(ssh2t.*ssh2t))

x1 = ncread('/Users/jearly/Documents/Data/ForcedDissipativeQGTurbulence/lowResSSHspatial.nc', 'x');
myssh1 = ncread('/Users/jearly/Documents/Data/ForcedDissipativeQGTurbulence/lowResSSHspatial.nc', 'ssh');
x2 = ncread('/Users/jearly/Documents/Data/ForcedDissipativeQGTurbulence/highResSSHspatial.nc', 'x');
myssh2 = ncread('/Users/jearly/Documents/Data/ForcedDissipativeQGTurbulence/highResSSHspatial.nc', 'ssh');

figure
plot(x1,myssh1(end,:), 'r')
hold on
plot(x2,myssh2(end,:), 'b')