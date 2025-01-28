 % - Topic: Energetics

The horizontally-averaged depth-integrated energy computed in the spatial domain is defined as,

$$
\frac{1}{2 L_x L_y} \int_0^{Lx} \int_0^{Ly} \int_{-L_z}^0 \left( u^2 + v^2 + w^2 + N^2 \eta^2 \right) dz
$$

for a non-hydrostatic transform and

$$
\frac{1}{2 L_x L_y} \int_0^{Lx} \int_0^{Ly} \int_{-L_z}^0 \left( u^2 + v^2 + N^2 \eta^2 \right) dz
$$

for a hydrostatic transform.

In code,

```matlab
if self.isHydrostatic == 1
    [u,v,eta] = self.variableWithName('u','v','eta');
    energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
else
    [u,v,w,eta] = self.variableWithName('u','v','w','eta');
    energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
end
```
