//
//  Quasigeostrophy2D.m
//  GLOceanKit
//
//  Created by Jeffrey J. Early on 10/2/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "Quasigeostrophy2D.h"

#define g 9.81
#define R 6.371e6


@implementation Quasigeostrophy2D

- (Quasigeostrophy2D *) initWithDimensions: (NSArray *) dims depth: (GLFloat) h latitude: (GLFloat) lat0 equation: (GLEquation *) equation
{
	if ((self=[super init])) {
		self.dimensions = dims;
		_h = h;
		_latitude = lat0;
		
		_f0 = 2 * 7.2921E-5 * sin( self.latitude*M_PI/180. );
		_beta = 2 * 7.2921E-5 * cos( self.latitude*M_PI/180. ) / R;
		
		// This is the only choice of parameters to completely nondimensionalize
		// the QGPVE on the beta-plane. f-plane has an additional freedom.
		_L_QG = sqrt(g*self.h)/self.f0; // m
		_T_QG = 1/(self.beta*self.L_QG); // s
		_N_QG = self.h*(self.beta*self.L_QG*self.L_QG)/sqrt(g*self.h); // m
		
		self.equation = equation;
		self.ssh = [GLFunction functionOfRealTypeWithDimensions: self.dimensions forEquation:self.equation];
		[self.ssh zero];
	}
	
	return self;
}

- (Quasigeostrophy2D *) initWithFile: (NSURL *) fileURL resolutionDoubling: (BOOL) shouldDouble
{
	GLEquation *equation = [[GLEquation alloc] init];
	GLNetCDFFile *restartFile = [[GLNetCDFFile alloc] initWithURL: fileURL forEquation: equation];
	
	// Now let's do a bunch of sanity checks to see if we really can read form this file.
	
	NSArray *requiredAttributes = @[@"is-anti-aliased", @"r", @"alpha", @"f_zeta", @"forcing-fraction-width", @"forcing-fraction", @"uses-beta", @"latitude", @"equivalent-depth", @"L_QG", @"N_QG", @"uses-spectral-vanishing-viscosity"];
	for (NSString *attribute in requiredAttributes) {
		if ( !(restartFile.globalAttributes[attribute]) ) {
			[NSException raise:@"InvalidRestartFileException" format: @"The restart file does not contain the attribute: %@.", attribute];
		}
	}
	
	GLFunction *allSSH = [restartFile variableWithName: @"SSH"];
	if ( !allSSH ) {
		[NSException raise:@"InvalidRestartFileException" format: @"The restart file does not contain an SSH history."];
	}
	
	if ( allSSH.dimensions.count != 3 ) {
		[NSException raise:@"InvalidRestartFileException" format: @"The SSH history does not contain the proper number of dimensions."];
	}
	
	// Extract the last SSH record
	NSArray *dims = allSSH.dimensions;
	NSMutableArray *ranges = [NSMutableArray array];
	ranges[0] = [NSValue valueWithRange: NSMakeRange([dims[0] nPoints]-1, 1)];
	ranges[1] = [NSValue valueWithRange: NSMakeRange(0, [dims[1] nPoints])];
	ranges[2] = [NSValue valueWithRange: NSMakeRange(0, [dims[2] nPoints])];
	GLFunction *ssh = [allSSH variableFromIndexRange: ranges];
	[ssh solve];
	
	// Nondimensionalize the ssh
	GLFloat L_QG = [restartFile.globalAttributes[@"L_QG"] doubleValue];
	GLFloat N_QG = [restartFile.globalAttributes[@"N_QG"] doubleValue];
	ssh = [ssh scaleVariableBy: 1./N_QG withUnits: @"unitless" dimensionsBy: 1./L_QG units: @"unitless"];
	
	//
	GLDimension *xDim = ssh.dimensions[0]; xDim.name = @"x";
	GLDimension *yDim = ssh.dimensions[1]; yDim.name = @"y";
	
	GLFunction *allPhase = [restartFile variableWithName: @"force_phase"];
	if (allPhase.dimensions.count == 3) {
		dims = allPhase.dimensions;
		ranges = [NSMutableArray array];
		ranges[0] = [NSValue valueWithRange: NSMakeRange([dims[0] nPoints]-1, 1)];
		ranges[1] = [NSValue valueWithRange: NSMakeRange(0, [dims[1] nPoints])];
		ranges[2] = [NSValue valueWithRange: NSMakeRange(0, [dims[2] nPoints])];
		allPhase = [allPhase variableFromIndexRange: ranges];
	}
	[allPhase solve];
	
	if (shouldDouble)
	{
		GLDimension *doubleX = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:2*(xDim.nPoints) domainMin:xDim.domainMin length:xDim.domainLength];
		doubleX.name = @"x";
		GLDimension *doubleY = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints:2*(yDim.nPoints) domainMin:yDim.domainMin length:yDim.domainLength];
		doubleY.name = @"y";
		
		GLFunction *sshX2 = [ssh projectOntoDimensions: @[doubleX, doubleY] usingSpectralBasis: @[@(kGLExponentialBasis), @(kGLExponentialBasis)]];
		
		xDim = doubleX;
		yDim = doubleY;
		ssh = [sshX2 spatialDomain];
	}
	
	GLFloat h = [restartFile.globalAttributes[@"equivalent-depth"] doubleValue];
	GLFloat lat = [restartFile.globalAttributes[@"latitude"] doubleValue];
	
	if ((self=[self initWithDimensions: @[xDim, yDim] depth: h latitude: lat equation: equation])) {
		self.ssh = ssh;
		self.shouldUseBeta = [restartFile.globalAttributes[@"uses-beta"] boolValue];
		self.shouldUseSVV = [restartFile.globalAttributes[@"uses-spectral-vanishing-viscosity"] boolValue];
		self.shouldAntiAlias = [restartFile.globalAttributes[@"is-anti-aliased"] boolValue];
		
		self.forcingFraction = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
		self.thermalDampingFraction = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
		self.frictionalDampingFraction = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
		self.forcingWidth = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
		self.f_zeta = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
		self.forcingDecorrelationTime = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
	}
	
	return self;
}

@end
