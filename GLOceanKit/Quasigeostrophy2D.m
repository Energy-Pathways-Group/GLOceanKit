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

@interface Quasigeostrophy2D ()
@property(strong) GLLinearTransform *inverseLaplacianMinusOne;
@property(strong) GLLinearTransform *diffJacobianX;
@property(strong) GLLinearTransform *diffJacobianY;
@property(strong) GLLinearTransform *diffLinear;
@property(strong) GLLinearTransform *tracerDamp;
@property BOOL shouldAdvancePhases;
@end

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
		self.shouldUseBeta = NO;
		self.shouldUseSVV = YES;
		self.shouldAntiAlias = NO;
		
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
	
	// Grab the dimensions (which won't be named because we just made everything unitless)
	GLDimension *xDim = ssh.dimensions[0]; xDim.name = @"x";
	GLDimension *yDim = ssh.dimensions[1]; yDim.name = @"y";
	
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
	
	// Now extract the phase of the forcing from the last record
	GLFunction *phi0 = [restartFile variableWithName: @"force_phase"];
	if (phi0.dimensions.count == 3) {
		dims = phi0.dimensions;
		ranges = [NSMutableArray array];
		ranges[0] = [NSValue valueWithRange: NSMakeRange([dims[0] nPoints]-1, 1)];
		ranges[1] = [NSValue valueWithRange: NSMakeRange(0, [dims[1] nPoints])];
		ranges[2] = [NSValue valueWithRange: NSMakeRange(0, [dims[2] nPoints])];
		phi0 = [phi0 variableFromIndexRange: ranges];
	}
	[phi0 solve];
	// And project onto the current spatial coordinates (regardless of whether or not we doubled).
	phi0 = [phi0 projectOntoDimensions: @[xDim, yDim] usingSpectralBasis: @[@(kGLExponentialBasis), @(kGLExponentialBasis)]];
	[phi0 solve];
	
	GLFloat h = [restartFile.globalAttributes[@"equivalent-depth"] doubleValue];
	GLFloat lat = [restartFile.globalAttributes[@"latitude"] doubleValue];
	
	if ((self=[self initWithDimensions: @[xDim, yDim] depth: h latitude: lat equation: equation])) {
		self.ssh = ssh;
		self.phi = phi0;
		self.shouldUseBeta = [restartFile.globalAttributes[@"uses-beta"] boolValue];
		self.shouldUseSVV = [restartFile.globalAttributes[@"uses-spectral-vanishing-viscosity"] boolValue];
		self.shouldAntiAlias = [restartFile.globalAttributes[@"is-anti-aliased"] boolValue];
		
		if (shouldDouble) {
			self.forcingFraction = 2*([restartFile.globalAttributes[@"forcing-fraction"] doubleValue]);
		} else {
			self.forcingFraction = [restartFile.globalAttributes[@"forcing-fraction"] doubleValue];
		}
		self.forcingWidth = [restartFile.globalAttributes[@"forcing-fraction-width"] doubleValue];
		self.f_zeta = [restartFile.globalAttributes[@"f_zeta"] doubleValue];
		self.forcingDecorrelationTime = [restartFile.globalAttributes[@"forcing-decorrelation-time"] doubleValue];
		self.thermalDampingFraction = [restartFile.globalAttributes[@"thermal-damping-fraction"] doubleValue];
		self.frictionalDampingFraction = [restartFile.globalAttributes[@"frictional-damping-fraction"] doubleValue];

	}
	
	return self;
}

-(void) setShouldAntiAlias:(BOOL)shouldAntiAlias {
	[GLBasisTransformOperation setShouldAntialias: shouldAntiAlias];
	_shouldAntiAlias = shouldAntiAlias;
}

- (void) createDifferentialOperators
{
	NSArray *spectralDimensions = [self.ssh dimensionsTransformedToBasis: self.ssh.differentiationBasis];
	
	GLLinearTransform *laplacian = [GLLinearTransform harmonicOperatorFromDimensions: spectralDimensions forEquation: self.equation];
	GLLinearTransform *laplacianMinusOne = [laplacian plus: @(-1.0)];
	self.inverseLaplacianMinusOne = [laplacianMinusOne inverse];
	
	GLLinearTransform *diff_xxx = [GLLinearTransform differentialOperatorWithDerivatives:@[@(3),@(0)] fromDimensions:spectralDimensions forEquation:self.equation];
	GLLinearTransform *diff_xyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(2)] fromDimensions:spectralDimensions forEquation:self.equation];
	GLLinearTransform *diff_xxy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(2),@(1)] fromDimensions:spectralDimensions forEquation:self.equation];
	GLLinearTransform *diff_yyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(0),@(3)] fromDimensions:spectralDimensions forEquation:self.equation];
	
	self.diffJacobianX = [diff_xxx plus: diff_xyy];
	self.diffJacobianY = [diff_xxy plus: diff_yyy];
	
	// I'm looking at Scott & Polvani, 2007 for this.
	NSArray *wavenumberDimensions = laplacianMinusOne.toDimensions;
	GLFunction *bigK;
	for (GLDimension *dim in wavenumberDimensions) {
		GLFunction *k = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: wavenumberDimensions forEquation:self.equation];
		bigK = !bigK ? [k multiply: k] : [bigK plus: [k multiply: k]];
	}
	bigK = [bigK sqrt];
	
	GLFloat k_max = [wavenumberDimensions[0] domainMin]+[wavenumberDimensions[0] domainLength];
	if ([GLBasisTransformOperation shouldAntialias]) {
		k_max = 2.0*k_max/3.0;
	}
	
	GLFloat k_f = k_max/self.forcingFraction;
	GLFloat k_width = [wavenumberDimensions[0] sampleInterval]*self.forcingWidth;
	
	GLFloat real_energy = pow(self.f_zeta,1.5)/(k_f*k_f);
	
	GLFloat k_alpha = [wavenumberDimensions[0] sampleInterval]*self.thermalDampingFraction;
	GLFloat alpha = pow(pow(self.f_zeta, 1.5)*pow(k_alpha,8.0)/pow(k_f, 2.0),0.3333);
	
	// Note that I'm adjusting the friction by what *really* happens. Our estimate is off, for some reason.
	GLFloat k_r = [wavenumberDimensions[0] sampleInterval]*self.frictionalDampingFraction;
	GLFloat r = 0.04*pow( pow(k_r, 2.0) * real_energy, 0.3333);
	
	// compute the svv cutoff, k_c
	GLFloat deltaK = [wavenumberDimensions[0] sampleInterval];
	GLFloat k_c = deltaK * pow(k_max/deltaK, 0.75);
	
	GLFloat u_rms = k_alpha > k_r ? pow(real_energy/k_alpha,0.333) : pow(real_energy/k_r,0.333);
	GLFloat nu = [self.dimensions[0] sampleInterval]*u_rms/2.0; // sqrt(2.) is very stable, 2. also seems to work.
	
	GLFloat k_nu = k_c;
	GLFloat k_min = k_c;
	GLFloat min = fabs(nu*exp( -pow((k_min-k_max)/(k_min-k_c),2.0) ) - sqrt(self.f_zeta)/(4*M_PI*M_PI*k_min*k_min));
	for (k_min=k_c; k_min<k_max; k_min += [wavenumberDimensions[0] sampleInterval]) {
		GLFloat a = fabs(nu*exp( -pow((k_min-k_max)/(k_min-k_c),2.0) ) - sqrt(self.f_zeta)/(4*M_PI*M_PI*k_min*k_min));
		if (a<min) {
			min=a;
			k_nu=k_min;
		}
	}
	
	// we aren't properly KVC here.
	_k_f = k_f;
	_k_nu = k_nu;
	_k_alpha = k_alpha;
	_k_r = k_r;
	_k_width = k_width;
	_k_max = k_max;
	_nu = nu;
	
	NSLog(@"k_nu=%f",k_nu);
	
	// The variable 'forcing' contains the magnitude of the forcing term for each wavenumber, but no phase information.
	GLVariable *forcing = [GLFunction functionOfRealTypeWithDimensions: wavenumberDimensions forEquation: self.equation];
	GLVariable *phaseSpeed;
	if (!self.phi) {
		// The random number generator will enforce the hermitian symmetry.
		self.phi = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: wavenumberDimensions forEquation: self.equation];
	}
	
	if (self.shouldForce)
	{
		// Now evenly distribute the intended force across all wavenumbers in the band
		GLFloat *f = [forcing pointerValue];
		GLFloat *kk = [bigK pointerValue];
		GLFloat totalWavenumbersInBand = 0;
		for (NSUInteger i=0; i<forcing.nDataPoints; i++) {
			if ( fabs(kk[i]-k_f) < k_width/2.0 ) totalWavenumbersInBand += 1.0;
		}
		// We are ignoring the negative frequencies for l. So we need to add those in as well.
		// This is not exactly double, but close enough.
		for (NSUInteger i=0; i<forcing.nDataPoints; i++) {
			if ( fabs(kk[i]-k_f) < k_width/2.0 ) f[i] = self.f_zeta/sqrt(2*totalWavenumbersInBand);
			else f[i] = 0.0;
		}
		
		if (self.forcingDecorrelationTime == HUGE_VAL) {
			// Infinite decorrelation time means the phases don't change.
			self.shouldAdvancePhases = NO;
		} else if (self.forcingDecorrelationTime > 0.0) {
			// Let the phases evolve with a give speed, but only if we chose a nonzero decorrelation time.
			phaseSpeed = [bigK scalarMultiply: 2*M_PI*self.T_QG/(k_f*self.forcingDecorrelationTime)];
			f = phaseSpeed.pointerValue;
			for (NSUInteger i=0; i<forcing.nDataPoints; i++) {
				if ( fabs(kk[i]-k_f) >= k_width/2.0 ) f[i] = 0.0;
			}
			self.shouldAdvancePhases = YES;
		} else {
			// In this case we will choose random values at each time step, so the phase change may as well be zero.
			self.shouldAdvancePhases = NO;
		}
		
	} else {
		[forcing zero];
		[self.phi zero];
	}
	
	/************************************************************************************************/
	/*		Create and cache the differential operators we will need								*/
	/************************************************************************************************/
	
	GLLinearTransform *harmonic = [GLLinearTransform harmonicOperatorOfOrder: 1 fromDimensions: spectralDimensions forEquation: self.equation];
	GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: spectralDimensions forEquation: self.equation];
	
	// first create the small scale damping operator, used for numerical stability.
	GLLinearTransform *numericalDamp;
	if (self.shouldUseSVV) {
		GLLinearTransform *svv = [GLLinearTransform spectralVanishingViscosityFilterWithDimensions: spectralDimensions scaledForAntialiasing: self.shouldAntiAlias forEquation: self.equation];
		numericalDamp = [[biharmonic times: @(nu)] times: svv];
		
		self.tracerDamp = [[harmonic times: @(nu)] times: svv];
	} else {
		GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: spectralDimensions forEquation: self.equation];
		numericalDamp = [biharmonic times: @(nu)];
		
		self.tracerDamp = [harmonic times: @(nu)];
	}
	
	// Now add the two large scale damping components.
	self.diffLinear = [[numericalDamp plus: @(alpha)] plus: [harmonic times: @(-r)]];
	
	if (self.shouldUseBeta) {
		GLLinearTransform *diff_x = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(0)] fromDimensions:spectralDimensions forEquation:self.equation];
		self.diffLinear = [self.diffLinear minus: diff_x];
	}
	
	/************************************************************************************************/
	/*		Estimate the time step																	*/
	/************************************************************************************************/
	
	GLDimension *xDim = self.dimensions[0];
	CGFloat cfl = 0.25;
	
	// Rounds dt to a number that evenly divides one day.
	GLFloat viscous_dt = cfl*xDim.sampleInterval * xDim.sampleInterval / (nu);
	GLFloat maxU = pow(alpha,-0.125)*pow(real_energy,3./8.);
	maxU = u_rms;
	GLFloat other_dt = cfl*xDim.sampleInterval/maxU;
	GLFloat dt = 1/(self.T_QG*ceil(1/(self.T_QG*other_dt)));
	
	if (shouldRestart) {
		GLVariable *v = [[ssh x] spatialDomain];
		GLVariable *u = [[ssh y] spatialDomain];
		GLVariable *speed = [[u times: u] plus: [v times: v]];
		[equation solveForVariable: speed];
		
		
		GLFloat U = sqrt([speed maxNow]);
		dt = cfl * xDim.sampleInterval / U;
	}
	
	NSLog(@"Reynolds number: %f", maxU*xDim.domainLength/nu);
	NSLog(@"v_dt: %f, other_dt: %f", viscous_dt*T_QG, other_dt*T_QG);
}

- (void) createIntegrationOperation
{

	
	NSMutableArray *yin = [NSMutableArray arrayWithObject: zeta];
	NSMutableArray *absTolerances = [NSMutableArray arrayWithObject: @(1e-6)];
	if (self.shouldAdvancePhases) {
		[yin addObject: self.phi];
		[absTolerances addObject: @(1e0)];
	}
	if (self.shouldAdvectFloats) {
		[yin addObjectsFromArray:@[self.xPosition, self.yPosition]];
		[absTolerances addObjectsFromArray:@[@(1e-3), @(1e-3)]];
	}
	if (self.shouldAdvectTracer) {
		for (GLFunction *tracer in self.tracers) {
			[yin addObject: [tracer frequencyDomain]];
			[absTolerances addObject: @(1e-3)];
		}
	}
	
	GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: yin stepSize: dt fFromTY:^(GLVariable *time, NSArray *yNew) {
		//		GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: yin stepSize: other_dt fFromTY:^(GLVariable *time, NSArray *yNew){
		NSUInteger iInput = 0;
		GLVariable *eta = [yNew[0] diff: @"inverseLaplacianMinusOne"];
		GLVariable *f = [[eta diff:@"diffLin"] plus: [[[[eta y] times: [eta diff: @"diffJacobianX"]] minus: [[eta x] times: [eta diff: @"diffJacobianY"]]] frequencyDomain]];
		
		if (shouldForce) {
			GLVariable *phase;
			if (forcingDecorrelationTime == HUGE_VAL) {
				phase = initialPhase;
			} else if (forcingDecorrelationTime == 0.0) {
				phase = [GLVariable variableWithRandomValuesBetween: -M_PI and: M_PI withDimensions: wavenumberDimensions forEquation: equation];
			} else {
				phase = yNew[++iInput];
			}
			f = [f plus: [forcing multiply: [[phase swapComplex] exponentiate]]];
		}
		
		NSMutableArray *fout = [NSMutableArray arrayWithObject: f];
		
		if (shouldAdvancePhases) {
			GLVariable *randomPhases = [GLVariable variableWithRandomValuesBetween: -M_PI and: M_PI withDimensions: wavenumberDimensions forEquation: equation];
			GLVariable *omega = [phaseSpeed multiply: [randomPhases dividedBy: [randomPhases abs]]];
			[fout addObject: omega];
		}
		
		if (shouldAdvectFloats) {
			NSArray *uv = @[[[[eta y] spatialDomain] negate], [[eta x] spatialDomain] ];
			NSArray *xy = @[yNew[++iInput], yNew[++iInput]];
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: uv secondOperand: xy];
			
			[fout addObjectsFromArray:@[interp.result[0], interp.result[1]]];
		}
		
		if (shouldAdvectTracer) {
			GLVariable *tracer = yNew[++iInput];
			GLVariable *fTracer = [[[[[eta y] times: [tracer x]] minus:[[eta x] times: [tracer y]] ] frequencyDomain] plus: [tracer diff: @"damp"]];
			
			[fout addObject: fTracer];
		}
		
		return fout;
		
	}];
	if ([[integrator class] isSubclassOfClass: [GLAdaptiveRungeKuttaOperation class]]) {
		[ (GLAdaptiveRungeKuttaOperation *) integrator setAbsoluteTolerance: absTolerances ];
	}
}

@end
