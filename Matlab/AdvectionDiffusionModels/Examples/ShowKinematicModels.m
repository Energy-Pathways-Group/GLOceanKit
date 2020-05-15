shouldSaveImages = 1;

eddy = TranslatingGaussian();
figure
eddy.plotStreamfunction(), hold on
eddy.plotVelocityField();

if shouldSaveImages == 1
    print('figures/kinematic_model_eddy.png','-dpng')
end

jet = MeanderingJet();
figure
jet.plotStreamfunction(), hold on
jet.plotVelocityField();

if shouldSaveImages == 1
    print('figures/kinematic_model_jet.png','-dpng')
end

cylinder = CylinderFlow();
figure
cylinder.plotStreamfunction(), hold on
cylinder.plotVelocityField();

if shouldSaveImages == 1
    print('figures/kinematic_model_cylinder.png','-dpng')
end

strain = LinearVelocityField(1e-6,0,0);
figure
strain.plotStreamfunction(), hold on
strain.plotVelocityField();

if shouldSaveImages == 1
    print('figures/kinematic_model_strain.png','-dpng')
end

strain = LinearVelocityField(0,0,1e-6);
figure
strain.plotStreamfunction(), hold on
strain.plotVelocityField();

if shouldSaveImages == 1
    print('figures/kinematic_model_vorticity.png','-dpng')
end