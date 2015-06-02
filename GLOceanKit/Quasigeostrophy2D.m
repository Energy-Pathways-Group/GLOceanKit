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

@property(strong) NSArray *dimensions;
@property(strong) NSArray *wavenumberDimensions;
@property(strong) GLMutableDimension *tDim;

@property(strong) GLEquation *equation;
@property(strong) GLRungeKuttaOperation *integrator;

@property(strong) GLLinearTransform *laplacian;
@property(strong) GLLinearTransform *laplacianMinusOne;
@property(strong) GLLinearTransform *diffJacobianX;
@property(strong) GLLinearTransform *diffJacobianY;
@property(strong) GLLinearTransform *diffLinear;
@property(strong) GLLinearTransform *tracerDamp;
@property BOOL shouldAdvancePhases;
@property BOOL isRestart;
@property GLFloat dt;

@property(strong) GLMutableVariable *sshHistory;
@property(strong) GLFunction *dimensionalForceMag;
@property(strong) GLMutableVariable *phaseHistory;
@property(strong) GLMutableVariable *sshFDHistory;
@property(strong) GLMutableVariable *rvHistory;
@property(strong) GLMutableVariable *forceHistory;
@property(strong) GLMutableVariable *xPositionHistory;
@property(strong) GLMutableVariable *yPositionHistory;
@property(strong) NSArray *tracerHistories;

@property GLFloat k_f;
@property GLFloat k_nu;
@property GLFloat k_alpha;
@property GLFloat k_r;
@property GLFloat k_width;
@property GLFloat k_max;
@property GLFloat alpha;
@property GLFloat r;
@property GLFloat nu;

- (void) createDifferentialOperators;
- (GLNetCDFFile *) createNetCDFFileAtURL: (NSURL *) outputURL;
- (void) createIntegrationOperation;

@end

@implementation Quasigeostrophy2D

- (Quasigeostrophy2D *) initWithDimensions: (NSArray *) dims depth: (GLFloat) h latitude: (GLFloat) lat0 equation: (GLEquation *) equation
{
	if ((self=[super init])) {
		
		if (dims.count != 2) {
			[NSException raise: @"InvalidDimensions" format: @"You must initialize with exactly two dimensions"];
		}
				
		_h = h;
		_latitude = lat0;
		
		_f0 = 2 * 7.2921E-5 * sin( self.latitude*M_PI/180. );
		_beta = 2 * 7.2921E-5 * cos( self.latitude*M_PI/180. ) / R;
		
		// This is the only choice of parameters to completely nondimensionalize
		// the QGPVE on the beta-plane. f-plane has an additional freedom.
		_L_QG = sqrt(g*self.h)/self.f0; // m
		_T_QG = 1/(self.beta*self.L_QG); // s
		_N_QG = self.h*(self.beta*self.L_QG*self.L_QG)/sqrt(g*self.h); // m
		
		GLDimension *xDim = [dims[0] scaledBy: 1/_L_QG translatedBy: 0.0 withUnits: @"unitless"];
		GLDimension *yDim = [dims[1] scaledBy: 1/_L_QG translatedBy: 0.0 withUnits: @"unitless"];
		self.dimensions = @[xDim, yDim];
		
		self.equation = equation;
		self.ssh = [GLFunction functionOfRealTypeWithDimensions: self.dimensions forEquation:self.equation];
		[self.ssh zero];
		self.shouldUseBeta = NO;
		self.shouldUseSVV = YES;
		self.shouldAntiAlias = NO;
		
	}
	
	return self;
}

- (Quasigeostrophy2D *) initWithFile: (NSURL *) fileURL resolutionDoubling: (BOOL) shouldDouble equation: (GLEquation *) equation
{
	GLNetCDFFile *restartFile = [[GLNetCDFFile alloc] initWithURL: fileURL forEquation: equation];
	
	// Now let's do a bunch of sanity checks to see if we really can read form this file.
	
	NSArray *requiredAttributes = @[@"is-anti-aliased", @"r", @"alpha", @"f_zeta", @"forcing-fraction-width", @"forcing-fraction", @"uses-beta", @"forcing-decorrelation-time", @"latitude", @"equivalent-depth", @"uses-spectral-vanishing-viscosity"];
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
	GLFloat L_QG = [restartFile.globalAttributes[@"length_scale"] doubleValue];
	GLFloat N_QG = [restartFile.globalAttributes[@"height_scale"] doubleValue];
    GLFloat T_QG = [restartFile.globalAttributes[@"time_scale"] doubleValue];
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
	
    NSLog(@"%g, %g", xDim.sampleInterval, [ssh.dimensions[0] sampleInterval ]);
    
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
    
    // Extract the magnitude of the force
    GLFunction *dimensionalForceMag = [restartFile variableWithName: @"force_magnitude"];
    GLFunction *forcing = [dimensionalForceMag scaleVariableBy: pow((T_QG),2.) withUnits: @"unitless" dimensionsBy: L_QG units: @"unitless"];
    forcing = [forcing projectOntoDimensions: @[xDim, yDim] usingSpectralBasis: @[@(kGLExponentialBasis), @(kGLExponentialBasis)]];
    [forcing solve];
    
	GLFloat h = [restartFile.globalAttributes[@"equivalent-depth"] doubleValue];
	GLFloat lat = [restartFile.globalAttributes[@"latitude"] doubleValue];
	
    GLDimension *dimXdim = [xDim scaledBy: L_QG translatedBy:0.0 withUnits:@"m"];
    GLDimension *dimYdim = [yDim scaledBy: L_QG translatedBy:0.0 withUnits:@"m"];
	if ((self=[self initWithDimensions: @[dimXdim, dimYdim] depth: h latitude: lat equation: equation])) {
		self.ssh = ssh;
		self.phi = phi0;
        self.forcing = forcing;
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
		self.isRestart = YES;
	}
	
	return self;
}

-(void) setShouldAntiAlias:(BOOL)shouldAntiAlias {
	[GLBasisTransformOperation setShouldAntialias: shouldAntiAlias];
	_shouldAntiAlias = shouldAntiAlias;
}

- (void) setShouldAdvectFloats:(BOOL)shouldAdvectFloats {
    // We are told to advect some floats, but the positions aren't set, so lets set the positions.
    if (shouldAdvectFloats && (!self.xPosition || !self.yPosition)) {
        GLDimension *xDim = self.dimensions[0];
        GLDimension *yDim = self.dimensions[1];
        GLDimension *xFloatDim = [[GLDimension alloc] initDimensionWithGrid: xDim.gridType nPoints: xDim.nPoints domainMin: xDim.domainMin length:xDim.domainLength];
        xFloatDim.name = @"x-float";
        GLDimension *yFloatDim = [[GLDimension alloc] initDimensionWithGrid: yDim.gridType nPoints: yDim.nPoints domainMin: yDim.domainMin length:yDim.domainLength];
        yFloatDim.name = @"y-float";
        
        NSArray *floatDims = @[xFloatDim, yFloatDim];
        self.xPosition = [GLFunction functionOfRealTypeFromDimension: xFloatDim withDimensions: floatDims forEquation: self.equation];
        self.yPosition = [GLFunction functionOfRealTypeFromDimension: yFloatDim withDimensions: floatDims forEquation: self.equation];
    }
    _shouldAdvectFloats = shouldAdvectFloats;
}

- (void) createDifferentialOperators
{
	self.wavenumberDimensions = [self.ssh dimensionsTransformedToBasis: self.ssh.differentiationBasis];
	
	self.laplacian = [GLLinearTransform harmonicOperatorFromDimensions: self.wavenumberDimensions forEquation: self.equation];
	self.laplacianMinusOne = [self.laplacian plus: @(-1.0)];
	self.inverseLaplacianMinusOne = [self.laplacianMinusOne inverse];
	
	GLLinearTransform *diff_xxx = [GLLinearTransform differentialOperatorWithDerivatives:@[@(3),@(0)] fromDimensions:self.wavenumberDimensions forEquation:self.equation];
	GLLinearTransform *diff_xyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(2)] fromDimensions:self.wavenumberDimensions forEquation:self.equation];
	GLLinearTransform *diff_xxy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(2),@(1)] fromDimensions:self.wavenumberDimensions forEquation:self.equation];
	GLLinearTransform *diff_yyy = [GLLinearTransform differentialOperatorWithDerivatives:@[@(0),@(3)] fromDimensions:self.wavenumberDimensions forEquation:self.equation];
	
	self.diffJacobianX = [diff_xxx plus: diff_xyy];
	self.diffJacobianY = [diff_xxy plus: diff_yyy];
	
	// I'm looking at Scott & Polvani, 2007 for this.
	GLFunction *bigK;
	for (GLDimension *dim in self.wavenumberDimensions) {
		GLFunction *k = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: self.wavenumberDimensions forEquation:self.equation];
		bigK = !bigK ? [k multiply: k] : [bigK plus: [k multiply: k]];
	}
	bigK = [bigK sqrt];
	
	GLFloat k_max = [self.wavenumberDimensions[0] domainMin]+[self.wavenumberDimensions[0] domainLength];
	if ([GLBasisTransformOperation shouldAntialias]) {
		k_max = 2.0*k_max/3.0;
	}
    
    if (self.f_zeta == 0) {
        NSLog(@"Warning: f_zeta is set to zero.");
    }
	
	self.k_f = k_max/self.forcingFraction;
	self.k_width = [self.wavenumberDimensions[0] sampleInterval]*self.forcingWidth;
	
	GLFloat real_energy = pow(self.f_zeta,1.5)/(self.k_f*self.k_f);
	
	self.k_alpha = [self.wavenumberDimensions[0] sampleInterval]*self.thermalDampingFraction;
	self.alpha = pow(pow(self.f_zeta, 1.5)*pow(self.k_alpha,8.0)/pow(self.k_f, 2.0),0.3333);
	
	// Note that I'm adjusting the friction by what *really* happens. Our estimate is off, for some reason.
	self.k_r = [self.wavenumberDimensions[0] sampleInterval]*self.frictionalDampingFraction;
	self.r = 0.04*pow( pow(self.k_r, 2.0) * real_energy, 0.3333);
	
	// compute the svv cutoff, k_c
	GLFloat deltaK = [self.wavenumberDimensions[0] sampleInterval];
	GLFloat k_c = deltaK * pow(k_max/deltaK, 0.75);
	
	GLFloat u_rms = self.k_alpha > self.k_r ? pow(real_energy/self.k_alpha,0.333) : pow(real_energy/self.k_r,0.333);
    if (u_rms==0.0 || isnan(u_rms)) {
        GLFunction *v = [[self.ssh x] spatialDomain];
        GLFunction *u = [[self.ssh y] spatialDomain];
        GLFunction *speed = [[u times: u] plus: [v times: v]];
        
        u_rms = sqrt([speed maxNow]);
    }
	self.nu = [self.dimensions[0] sampleInterval]*u_rms; // sqrt(2.) is very stable, 2. also seems to work.
	
    NSLog(@"sample interval: %g, u_rms: %g, real_energy: %g",[self.dimensions[0] sampleInterval],u_rms ,real_energy);
    
	self.k_nu = k_c;
	GLFloat k_min = k_c;
	GLFloat min = fabs(self.nu*exp( -pow((k_min-k_max)/(k_min-k_c),2.0) ) - sqrt(self.f_zeta)/(4*M_PI*M_PI*k_min*k_min));
	for (k_min=k_c; k_min<k_max; k_min += [self.wavenumberDimensions[0] sampleInterval]) {
		GLFloat a = fabs(self.nu*exp( -pow((k_min-k_max)/(k_min-k_c),2.0) ) - sqrt(self.f_zeta)/(4*M_PI*M_PI*k_min*k_min));
		if (a<min) {
			min=a;
			self.k_nu=k_min;
		}
	}
	

	self.k_max = k_max;
		
    [self buildForcingFunction];
	
	/************************************************************************************************/
	/*		Create and cache the differential operators we will need								*/
	/************************************************************************************************/
	
	GLLinearTransform *harmonic = [GLLinearTransform harmonicOperatorOfOrder: 1 fromDimensions: self.wavenumberDimensions forEquation: self.equation];
	GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: self.wavenumberDimensions forEquation: self.equation];
	
	// first create the small scale damping operator, used for numerical stability.
	GLLinearTransform *numericalDamp;
	if (self.shouldUseSVV) {
		GLLinearTransform *svv = [GLLinearTransform spectralVanishingViscosityFilterWithDimensions: self.wavenumberDimensions scaledForAntialiasing: self.shouldAntiAlias forEquation: self.equation];
		numericalDamp = [[biharmonic times: @(self.nu)] times: svv];
		
		self.tracerDamp = [[harmonic times: @(self.nu)] times: svv];
	} else {
		GLLinearTransform *biharmonic = [GLLinearTransform harmonicOperatorOfOrder: 2 fromDimensions: self.wavenumberDimensions forEquation: self.equation];
		numericalDamp = [biharmonic times: @(self.nu)];
		
		self.tracerDamp = [harmonic times: @(self.nu)];
	}
	
	// Now add the two large scale damping components.
	self.diffLinear = [[numericalDamp plus: @(self.alpha)] plus: [harmonic times: @(-self.r)]];
	
	if (self.shouldUseBeta) {
		GLLinearTransform *diff_x = [GLLinearTransform differentialOperatorWithDerivatives:@[@(1),@(0)] fromDimensions:self.wavenumberDimensions forEquation:self.equation];
		self.diffLinear = [self.diffLinear minus: diff_x];
	}
	
	/************************************************************************************************/
	/*		Estimate the time step																	*/
	/************************************************************************************************/
	
	GLDimension *xDim = self.dimensions[0];
	CGFloat cfl = 0.25;
	
	// Rounds dt to a number that evenly divides one day.
	GLFloat viscous_dt = cfl*xDim.sampleInterval * xDim.sampleInterval / (self.nu);
	GLFloat other_dt = cfl*xDim.sampleInterval/u_rms;
	//self.dt = 1/(self.T_QG*ceil(1/(self.T_QG*other_dt)));
#warning come back and fix the time step issue.
    self.dt = other_dt;
	
    if (isnan(self.dt)) {
		GLFunction *v = [[self.ssh x] spatialDomain];
		GLFunction *u = [[self.ssh y] spatialDomain];
		GLFunction *speed = [[u times: u] plus: [v times: v]];
		
		GLFloat U = sqrt([speed maxNow]);
        if (U==0.0) U=1.0;
		self.dt = cfl * xDim.sampleInterval / U;
	}
	
	NSLog(@"Reynolds number: %f", u_rms*xDim.domainLength/self.nu);
	NSLog(@"v_dt: %f, other_dt: %f", viscous_dt*self.T_QG, other_dt*self.T_QG);
    NSLog(@"dt: %f", self.dt);
}

- (void) buildForcingFunction
{
    if (self.shouldForce) {
        
        if (!self.phi) {
            // The random number generator will enforce the hermitian symmetry.
            self.phi = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: self.wavenumberDimensions forEquation: self.equation];
            // This line is pivotal: we must remove the dependency to the random number generator operation before using this.
            [self.phi solve];
        }
        
        // The variable 'forcing' contains the magnitude of the forcing term for each wavenumber, but no phase information.
        if (!self.forcing) {
            self.forcing = [GLFunction functionOfRealTypeWithDimensions: self.wavenumberDimensions forEquation: self.equation];
            
            GLFunction *bigK;
            for (GLDimension *dim in self.wavenumberDimensions) {
                GLFunction *k = [GLFunction functionOfRealTypeFromDimension: dim withDimensions: self.wavenumberDimensions forEquation:self.equation];
                bigK = !bigK ? [k multiply: k] : [bigK plus: [k multiply: k]];
            }
            bigK = [bigK sqrt];
            
            // Now evenly distribute the intended force across all wavenumbers in the band
            GLFloat *f = [self.forcing pointerValue];
            GLFloat *kk = [bigK pointerValue];
            GLFloat totalWavenumbersInBand = 0;
            for (NSUInteger i=0; i<self.forcing.nDataPoints; i++) {
                if ( fabs(kk[i]-self.k_f) < self.k_width/2.0 ) totalWavenumbersInBand += 1.0;
            }
            // We are ignoring the negative frequencies for l. So we need to add those in as well.
            // This is not exactly double, but close enough.
            for (NSUInteger i=0; i<self.forcing.nDataPoints; i++) {
                if ( fabs(kk[i]-self.k_f) < self.k_width/2.0 ) f[i] = self.f_zeta/sqrt(2*totalWavenumbersInBand);
                else f[i] = 0.0;
            }
            
            if (self.forcingDecorrelationTime == HUGE_VAL) {
                // Infinite decorrelation time means the phases don't change.
                self.shouldAdvancePhases = NO;
            } else if (self.forcingDecorrelationTime > 0.0) {
                // Let the phases evolve with a give speed, but only if we chose a nonzero decorrelation time.
                self.phaseSpeed = [bigK scalarMultiply: 2*M_PI*self.T_QG/(self.k_f*self.forcingDecorrelationTime)];
                f = self.phaseSpeed.pointerValue;
                for (NSUInteger i=0; i<self.forcing.nDataPoints; i++) {
                    if ( fabs(kk[i]-self.k_f) >= self.k_width/2.0 ) f[i] = 0.0;
                }
                self.shouldAdvancePhases = YES;
            } else {
                // In this case we will choose random values at each time step, so the phase change may as well be zero.
                self.shouldAdvancePhases = NO;
            }
        }
    } else {
        self.forcing = nil;
        self.phi = nil;
    }
}

- (void) createIntegrationOperation
{
	GLFunction *zeta = [self.ssh differentiateWithOperator: self.laplacianMinusOne];
    GLFunction *zetaSD = [zeta spatialDomain];
	
	NSMutableArray *yin = [NSMutableArray arrayWithObject: zeta];
	NSMutableArray *absTolerances = [NSMutableArray arrayWithObject: @(1e-1)];
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
	
    __weak Quasigeostrophy2D *qg = self;
	self.fBlock = ^(GLVariable *time, NSArray *yNew) {
		NSUInteger iInput = 0;
		GLFunction *eta = [yNew[0] differentiateWithOperator: qg.inverseLaplacianMinusOne];
		//GLFunction *f = [[eta differentiateWithOperator: qg.diffLinear] plus: [[[[eta y] times: [eta differentiateWithOperator: qg.diffJacobianX]] minus: [[eta x] times: [eta differentiateWithOperator: qg.diffJacobianY]]] frequencyDomain]];
        
        GLFunction *eta_y = [eta y];
        GLFunction *eta_nabla_x = [eta differentiateWithOperator: qg.diffJacobianX];
        GLFunction *eta_x = [eta x];
        GLFunction *eta_nabla_y = [eta differentiateWithOperator: qg.diffJacobianY];
        
        GLFunction *eta_y_sd = [eta_y spatialDomain];
        GLFunction *eta_nabla_x_sd = [eta_nabla_x spatialDomain];
        GLFunction *eta_x_sd = [eta_x spatialDomain];
        GLFunction *eta_nabla_y_sd = [eta_nabla_y spatialDomain];
        
        GLFunction *a = [eta_y_sd times: eta_nabla_x_sd];
        GLFunction *b = [eta_x_sd times: eta_nabla_y_sd];
        
        GLFunction *f = [[eta differentiateWithOperator: qg.diffLinear] plus: [[a minus: b] frequencyDomain]];
        

		if (qg.shouldForce) {
			GLFunction *phase;
			if (qg.forcingDecorrelationTime == HUGE_VAL) {
				phase = qg.phi;
			} else if (qg.forcingDecorrelationTime == 0.0) {
				phase = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: qg.wavenumberDimensions forEquation: qg.equation];
			} else {
				phase = yNew[++iInput];
			}
			f = [f plus: [qg.forcing multiply: [[phase swapComplex] exponentiate]]];            
		}
		
		NSMutableArray *fout = [NSMutableArray arrayWithObject: f];
		
		if (qg.shouldAdvancePhases) {
			GLFunction *randomPhases = [GLFunction functionWithRandomValuesBetween: -M_PI and: M_PI withDimensions: qg.wavenumberDimensions forEquation: qg.equation];
			GLFunction *omega = [qg.phaseSpeed multiply: [randomPhases dividedBy: [randomPhases abs]]];
			[fout addObject: omega];
		}
		
		if (qg.shouldAdvectFloats) {
			NSArray *uv = @[[[[eta y] spatialDomain] negate], [[eta x] spatialDomain] ];
			NSArray *xy = @[yNew[++iInput], yNew[++iInput]];
			GLSimpleInterpolationOperation *interp = [[GLSimpleInterpolationOperation alloc] initWithFirstOperand: uv secondOperand: xy];
			
			[fout addObjectsFromArray:@[interp.result[0], interp.result[1]]];
		}
		
		if (qg.shouldAdvectTracer) {
			for (NSUInteger iTracer=0; iTracer<qg.tracers.count; iTracer++) {
				GLFunction *tracer = yNew[++iInput];
				GLFunction *fTracer = [[[[[eta y] times: [tracer x]] minus:[[eta x] times: [tracer y]] ] frequencyDomain] plus: [tracer differentiateWithOperator: qg.tracerDamp]];
				
				[fout addObject: fTracer];
			}
		}
		
		return fout;
		
	};
    GLAdaptiveRungeKuttaOperation *integrator = [GLAdaptiveRungeKuttaOperation rungeKutta23AdvanceY: yin stepSize: self.dt fFromTY: self.fBlock];
    //GLRungeKuttaOperation *integrator = [GLRungeKuttaOperation rungeKutta4AdvanceY: yin stepSize: self.dt fFromTY: self.fBlock];
    if ([[integrator class] isSubclassOfClass: [GLAdaptiveRungeKuttaOperation class]]) {
		[ (GLAdaptiveRungeKuttaOperation *) integrator setAbsoluteTolerance: absTolerances ];
	}
    
    self.integrator = integrator;
    self.yin = yin;
    self.absoluteTolerance = absTolerances;
}

- (void) addMetadataToNetCDFFile: (GLNetCDFFile *) netcdfFile;
{
	GLFloat wavenumberScale = 1./(self.L_QG);
	
	[netcdfFile setGlobalAttribute: @(self.shouldUseBeta) forKey:@"uses-beta"];
	[netcdfFile setGlobalAttribute: @(self.shouldAntiAlias) forKey:@"is-anti-aliased"];
	[netcdfFile setGlobalAttribute: @(self.shouldUseSVV) forKey:@"uses-spectral-vanishing-viscosity"];
	[netcdfFile setGlobalAttribute: @(self.forcingFraction) forKey:@"forcing-fraction"];
	[netcdfFile setGlobalAttribute: @(self.thermalDampingFraction) forKey:@"thermal-damping-fraction"];
	[netcdfFile setGlobalAttribute: @(self.frictionalDampingFraction) forKey:@"frictional-damping-fraction"];
	[netcdfFile setGlobalAttribute: @(self.forcingWidth) forKey:@"forcing-fraction-width"];
	[netcdfFile setGlobalAttribute: @(self.f_zeta) forKey:@"f_zeta"];
	[netcdfFile setGlobalAttribute: @(self.h) forKey:@"equivalent-depth"];
	[netcdfFile setGlobalAttribute: @(self.alpha) forKey:@"alpha"];
	[netcdfFile setGlobalAttribute: @(self.r) forKey:@"r"];
	[netcdfFile setGlobalAttribute: @(self.forcingDecorrelationTime) forKey:@"forcing-decorrelation-time"];
	[netcdfFile setGlobalAttribute: @(self.latitude) forKey:@"latitude"];
	
	// The rest of these can be derived from the others, but it's convinient to save them anyway.
	[netcdfFile setGlobalAttribute: @(self.k_f*wavenumberScale) forKey:@"forcing_wavenumber"];
	[netcdfFile setGlobalAttribute: @(self.k_nu*wavenumberScale) forKey:@"viscous_wavenumber"];
	[netcdfFile setGlobalAttribute: @(self.k_alpha*wavenumberScale) forKey:@"thermal_damping_wavenumber"];
	[netcdfFile setGlobalAttribute: @(self.k_r*wavenumberScale) forKey:@"frictional_damping_wavenumber"];
	[netcdfFile setGlobalAttribute: @(self.k_width*wavenumberScale) forKey:@"forcing_width"];
	[netcdfFile setGlobalAttribute: @(self.k_max*wavenumberScale) forKey:@"max_resolved_wavenumber"];
	[netcdfFile setGlobalAttribute: @(self.nu*self.L_QG*self.L_QG/(self.T_QG)) forKey:@"nu"];
	[netcdfFile setGlobalAttribute: @(self.N_QG) forKey:@"height_scale"];
	[netcdfFile setGlobalAttribute: @(self.L_QG) forKey:@"length_scale"];
	[netcdfFile setGlobalAttribute: @(1/(self.T_QG)) forKey:@"vorticity_scale"];
	[netcdfFile setGlobalAttribute: @(self.T_QG) forKey:@"time_scale"];
}

- (GLNetCDFFile *) createNetCDFFileAtURL: (NSURL *) outputURL
{
	GLNetCDFFile *netcdfFile = [[GLNetCDFFile alloc] initWithURL: outputURL forEquation: self.equation overwriteExisting: YES];
	[self addMetadataToNetCDFFile: netcdfFile];
	
	self.tDim = [[GLMutableDimension alloc] initWithPoints: @[@(0.0)]];
	self.tDim.name = @"time";
	self.tDim.units = @"s";
	
	GLFloat wavenumberScale = 1./(self.L_QG);
	
	GLFunction *dimensionalSSH = [[self.ssh spatialDomain] scaleVariableBy: self.N_QG withUnits: @"m" dimensionsBy: self.L_QG units: @"m"];
	self.sshHistory = [dimensionalSSH variableByAddingDimension: self.tDim];
	self.sshHistory.name = @"SSH";
	self.sshHistory = [netcdfFile addVariable: self.sshHistory];
	
    if (self.forcing) {
        GLFunction *dimensionalForceMag = [self.forcing scaleVariableBy: pow((self.T_QG),-2.) withUnits: @"1/s^2" dimensionsBy: wavenumberScale units: @"cycles/m"];
        dimensionalForceMag.name = @"force_magnitude";
        self.dimensionalForceMag = [netcdfFile addVariable: dimensionalForceMag];
    }
    
    if (self.phi) {
        GLFunction *dimensionalPhase = [self.phi scaleVariableBy: 1.0 withUnits: @"radians" dimensionsBy: wavenumberScale units: @"cycles/m"];
        dimensionalPhase.name = @"force_phase";
        if (self.shouldAdvancePhases) {
            self.phaseHistory = [dimensionalPhase variableByAddingDimension: self.tDim];
            self.phaseHistory.name = @"force_phase";
            self.phaseHistory = [netcdfFile addVariable: self.phaseHistory];
        } else {
            [netcdfFile addVariable: dimensionalPhase];
        }
    }
	
	if (self.shouldWriteSSHFD) {
		GLFunction *dimensionalSSHFD = [[self.ssh frequencyDomain] scaleVariableBy: self.N_QG withUnits: @"m^3" dimensionsBy: wavenumberScale units: @"cycles/m"];
		self.sshFDHistory = [dimensionalSSHFD variableByAddingDimension: self.tDim];
		self.sshFDHistory.name = @"SSH_FD";
		self.sshFDHistory = [netcdfFile addVariable: self.sshFDHistory];
	}
	
	if (self.shouldWriteRV) {
		GLFloat rvScale = (g/self.f0)*(self.N_QG)/pow(self.L_QG,2.0);
		GLFunction *dimensionalRV = [[[self.ssh differentiateWithOperator: self.laplacian] spatialDomain] scaleVariableBy: rvScale withUnits: @"1/s" dimensionsBy: self.L_QG units: @"m"];
		self.rvHistory = [dimensionalRV variableByAddingDimension: self.tDim];
		self.rvHistory.name = @"RV";
		self.rvHistory = [netcdfFile addVariable: self.rvHistory];
	}
	
	if (self.shouldWriteForce) {
		GLFunction *dimensionalForce = [[[self.forcing multiply: [[self.phi swapComplex] exponentiate]] spatialDomain] scaleVariableBy: pow((self.T_QG),-2.) withUnits: @"1/s^2" dimensionsBy: self.L_QG units: @"m"];
		self.forceHistory = [dimensionalForce variableByAddingDimension: self.tDim];
		self.forceHistory.name = @"force";
		self.forceHistory = [netcdfFile addVariable: self.forceHistory];
	}
	
	if (self.shouldAdvectFloats) {
		GLFunction *dimensionalXPosition = [self.xPosition scaleVariableBy: self.L_QG withUnits: @"m" dimensionsBy: self.L_QG units: @"m"];
		self.xPositionHistory = [dimensionalXPosition variableByAddingDimension: self.tDim];
		self.xPositionHistory.name = @"x-position";
		self.xPositionHistory = [netcdfFile addVariable: self.xPositionHistory];
		
		GLFunction *dimensionalYPosition = [self.yPosition scaleVariableBy: self.L_QG withUnits: @"m" dimensionsBy: self.L_QG units: @"m"];
		self.yPositionHistory = [dimensionalYPosition variableByAddingDimension: self.tDim];
		self.yPositionHistory.name = @"y-position";
		self.yPositionHistory = [netcdfFile addVariable: self.yPositionHistory];
	}
	
	if (self.shouldAdvectTracer) {
		NSMutableArray *tracerHistories = [NSMutableArray array];
		for (GLFunction *tracer in self.tracers) {
			GLMutableVariable *tracerHistory = [tracer variableByAddingDimension: self.tDim];
			tracerHistory.name = @"tracer";
			tracerHistory = [netcdfFile addVariable: tracerHistory];
			[tracerHistories addObject: tracerHistory];
		}
		self.tracerHistories = tracerHistories;
	}
	
	return netcdfFile;
}


- (void) runSimulationToTime: (GLFloat) maxTime
{
    [self createDifferentialOperators];
    [self createIntegrationOperation];
    
    GLNetCDFFile *netcdfFile;
    if (self.outputFile) {
        netcdfFile = [self createNetCDFFileAtURL: self.outputFile];
    }
    
    GLFloat wavenumberScale = 1./(self.L_QG);
    GLFloat rvScale = (g/self.f0)*(self.N_QG)/pow(self.L_QG,2.0);
    
	
	GLFloat outputIntervalND = self.outputInterval/self.T_QG;
	GLFloat maxTimeND = maxTime/self.T_QG;
	GLFloat stepND =self.outputInterval/self.T_QG;
	NSLog(@"Starting simulation...running to maxTime %f",maxTimeND);
    for (GLFloat time = self.outputInterval/self.T_QG; time < maxTime/self.T_QG; time += self.outputInterval/self.T_QG)
    {
        @autoreleasepool {
            NSUInteger iOut = 0;
            NSArray *yout = [self.integrator stepForwardToTime: time];
            NSLog(@"Logging day: %f, last step size: %f seconds.", (time*self.T_QG/86400.), self.integrator.lastStepSize*self.T_QG);
            
            [self.tDim addPoint: @(time*self.T_QG)];
            
            GLFunction *etaFD = [yout[0] differentiateWithOperator: self.inverseLaplacianMinusOne];
            self.ssh = [[etaFD spatialDomain] scaleVariableBy: self.N_QG withUnits: @"m" dimensionsBy: self.L_QG units: @"m"];
            
            [self.sshHistory concatenateWithLowerDimensionalVariable: self.ssh alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
            [netcdfFile waitUntilAllOperationsAreFinished];
            
            if (self.shouldWriteSSHFD) {
                GLFunction *etaFDDim = [etaFD scaleVariableBy: self.N_QG withUnits: @"m^3" dimensionsBy: wavenumberScale units: @"cycles/m"];
                [self.sshFDHistory concatenateWithLowerDimensionalVariable: etaFDDim alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
            
            if (self.shouldWriteRV) {
                GLFunction *rv2 = [[[self.laplacian transform: etaFD] spatialDomain] scaleVariableBy: rvScale withUnits: @"1/s" dimensionsBy: self.L_QG units: @"m"];
                [self.rvHistory concatenateWithLowerDimensionalVariable: rv2 alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
            
            if (self.shouldAdvancePhases) {
                GLFunction *phase = [yout[++iOut] scaleVariableBy: 1.0 withUnits: @"radians" dimensionsBy: wavenumberScale units: @"cycles/m"];
                [self.phaseHistory concatenateWithLowerDimensionalVariable: phase alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
                
                if (self.shouldWriteForce) {
                    GLFunction *dimensionalForce = [[[self.forcing multiply: [[yout[iOut] swapComplex] exponentiate]] spatialDomain] scaleVariableBy: pow((self.T_QG),-2.) withUnits: @"1/s^2" dimensionsBy: self.L_QG units: @"m"];
                    [self.forceHistory concatenateWithLowerDimensionalVariable: dimensionalForce alongDimensionAtIndex:0  toIndex:(self.tDim.nPoints-1)];
                }
            }
            
            if (self.shouldAdvectFloats) {
                self.xPosition = [yout[++iOut] scaleVariableBy: self.L_QG withUnits: @"m" dimensionsBy: self.L_QG units: @"m"];
                self.yPosition = [yout[++iOut] scaleVariableBy: self.L_QG withUnits: @"m" dimensionsBy: self.L_QG units: @"m"];
                [self.xPositionHistory concatenateWithLowerDimensionalVariable: self.xPosition alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
                [netcdfFile waitUntilAllOperationsAreFinished];
                [self.yPositionHistory concatenateWithLowerDimensionalVariable: self.yPosition alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
                [netcdfFile waitUntilAllOperationsAreFinished];
            }
            
            if (self.shouldAdvectTracer) {
                for (GLMutableVariable *tracerHistory in self.tracerHistories) {
                    GLFunction *tracer = [yout[++iOut] spatialDomain];
                    [tracerHistory concatenateWithLowerDimensionalVariable: tracer alongDimensionAtIndex:0 toIndex: (self.tDim.nPoints-1)];
                    [netcdfFile waitUntilAllOperationsAreFinished];
                }
            }
        }
    }
    
    NSLog(@"Close the NetCDF file and wrap up");
    
    [self.equation waitUntilAllOperationsAreFinished];
    
    // The NetCDF file may still be writing data. We need to make sure it finishes before we exit.
    [netcdfFile waitUntilAllOperationsAreFinished];
    [netcdfFile close];
}

@end
