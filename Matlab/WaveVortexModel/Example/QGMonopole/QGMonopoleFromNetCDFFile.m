model = WaveVortexModel.modelFromFile('QGMonopole.nc');

model.setupIntegrator(timeStepConstraint="advective", outputInterval=86400);
model.createNetCDFFileForModelOutput('QGMonopole-restart.nc',shouldOverwriteExisting=1);
model.integrateToTime(150*86400);

wvt = model.wvt;

figure, pcolor(wvt.x,wvt.y,wvt.ssh.'), shading interp