% program to divide cyprus data into separate files for each processor.
nprocs=16;
u=importdata('cyprus_um.txt');
v=importdata('cyprus_vm.txt');
ny=size(u,2)
nx=ny;
nz=size(u,1)/nx;
u=reshape(u,[nx nz ny]);
v=reshape(v,[nx nz ny]);
u=permute(u,[1 3 2]);
v=permute(v,[1 3 2]);
locnz=nz/nprocs;
u2d=zeros(nx,ny,locnz);
v2d=zeros(nx,ny,locnz);
for n=1:nprocs
myid=n-1
u2d=u(:,:,myid*locnz+1:(myid+1)*locnz);
v2d=v(:,:,myid*locnz+1:(myid+1)*locnz);
myid*locnz+1:(myid+1)*locnz
Flnu=sprintf('cyprus_u_%02d.txt',myid)
Flnv=sprintf('cyprus_v_%02d.txt',myid)
for k=1:locnz
u2=u2d(:,:,k);
v2=v2d(:,:,k);
save(Flnu,'u2','-ASCII','-DOUBLE','-APPEND')
save(Flnv,'v2','-ASCII','-DOUBLE','-APPEND')
end
end

