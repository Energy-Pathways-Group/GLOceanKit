The x coordinate is periodic, which means that
```matlab
dx = Lx/Nx;
x = dx*(0:Nx-1)';
```

Note that this means that it is NOT true that Lx=x(end)-x(1), but in fact you need an extra grid point, i.e., it IS true that Lx = dx + x(end)-x(1). This is the usual grid for Fourier Transforms.
