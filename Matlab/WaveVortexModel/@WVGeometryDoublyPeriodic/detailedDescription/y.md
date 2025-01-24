- Topic: Domain attributes — Grid — Spatial
- nav_order: 2

The values `Ly` and `Ny` are set during initialization from which the `y` coordinate is derived.

The y coordinate is periodic, which means that
```matlab
dy = Ly/Ny;
y = dy*(0:Ny-1)';
```

Note that this means that it is NOT true that Ly=y(end)-y(1), but in fact you need an extra grid point, i.e., it IS true that Ly = dy + y(end)-y(1). This is the usual grid for Fourier Transforms.
