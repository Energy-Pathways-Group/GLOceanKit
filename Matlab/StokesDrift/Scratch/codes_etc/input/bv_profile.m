% specify bv profiles to the north and south of domain
bv2max_south=0.8e-6; % (radians/s)^2
zthermo1_south = -2000.; % depth of first thermocline to the south
decay_up_south = 2000.; % decay rate above first thermocline in south
decay_down_south = 1000.; % decay rate below first thermocline in south

bv2max_north=0.8e-6; % (radians/s)^2
zthermo1_north = -2000.; % depth of first thermocline to the north
decay_up_north = 2000.; % decay rate above first thermocline in north
decay_down_north = 1000.; % decay rate below first thermocline in north

z1_south = -500.; % depth of first (shallow) bump in the south
z1_north = -750.; % depth of first (shallow) bump in the north
bv2max1_south = 2.5e-5; % bv^2 for shallow bump in south
bv2max1_north = 0.9e-5; % bv^2 for shallow bump in south
decay1_up_south  = 

sigma_bottom = 27.915 ; % bottom value of sigma

% calculating the squared Brunt-Vaisala frequency
% define the vertical coordinate
Lz = 4000.% domain depth [m]
nz = 512; %vertical resolution
dz = float(Lz/nz);
z=zeros(nz,1);

for k=0:nz
z(k)=k*dz;
if  z(k) > zthermo1_south  
arg = (z(k)-zthermo1_south)/decay_up_south^2 ;
elseif z(k) < zthermo1_south
arg = (z(k)-zthermo1_south)/decay_down_south^2 ;
end
bvsn






