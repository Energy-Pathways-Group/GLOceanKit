z = linspace(-5000, 0, 500);


imAnalytical = InternalModesExponentialStratification([5.2e-3 1300], [-5000 0], z, 33,'nModes',128);
imAnalytical.normalization = Normalization.kConstant;

profile on
for i=1:10
    [F_analytical,G_analytical,h_analytical] = imAnalytical.ModesAtFrequency( 0.2*imAnalytical.N0 );
end

profile viewer