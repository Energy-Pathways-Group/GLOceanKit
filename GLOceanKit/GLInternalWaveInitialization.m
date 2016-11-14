//
//  GLInternalWaveInitialization.m
//  InternalWaves
//
//  Created by Jeffrey J. Early on 1/14/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import "GLInternalWaveInitialization.h"
#import "GLInternalModes.h"

#define g 9.81

@interface GLInternalWaveInitialization ()

@property(strong) GLDimension *kDim;
@property(strong) GLDimension *lDim;
@property(strong) GLDimension *modeDim;
@property(strong) NSMutableArray *horizontalDimensions;

- (void) generateWavePhasesFromPositive: (GLFunction *) G_plus negative: (GLFunction *) G_minus;

@end

static NSString *GLInternalWaveMaximumModesKey = @"GLInternalWaveMaximumModesKey";
static NSString *GLInternalWaveEquationKey = @"GLInternalWaveEquationKey";
static NSString *GLInternalWaveInternalModeKey = @"GLInternalWaveInternalModeKey";
static NSString *GLInternalWaveFullDimensionsKey = @"GLInternalWaveFullDimensionsKey";
static NSString *GLInternalWaveSpectralDimensionsKey = @"GLInternalWaveSpectralDimensionsKey";
static NSString *GLInternalWaveLatitudeKey = @"GLInternalWaveLatitudeKey";
static NSString *GLInternalWaveF0Key = @"GLInternalWaveF0Key";
static NSString *GLInternalWaveRhoKey = @"GLInternalWaveRhoKey";
static NSString *GLInternalWaveN2Key = @"GLInternalWaveN2Key";
static NSString *GLInternalWaveEigenfrequenciesKey = @"GLInternalWaveEigenfrequenciesKey";
static NSString *GLInternalWaveEigendepthsKey = @"GLInternalWaveEigendepthsKey";
static NSString *GLInternalWaveRossbyRadiiKey = @"GLInternalWaveRossbyRadiiKey";
static NSString *GLInternalWaveSKey = @"GLInternalWaveSKey";
static NSString *GLInternalWaveSprimeKey = @"GLInternalWaveSprimeKey";
static NSString *GLInternalWaveZetaPlusKey = @"GLInternalWaveZetaPlusKey";
static NSString *GLInternalWaveZetaMinusKey = @"GLInternalWaveZetaMinusKey";
static NSString *GLInternalWaveRhoPlusKey = @"GLInternalWaveRhoPlusKey";
static NSString *GLInternalWaveRhoMinusKey = @"GLInternalWaveRhoMinusKey";
static NSString *GLInternalWaveUPlusKey = @"GLInternalWaveUPlusKey";
static NSString *GLInternalWaveUMinusKey = @"GLInternalWaveUMinusKey";
static NSString *GLInternalWaveVPlusKey = @"GLInternalWaveVPlusKey";
static NSString *GLInternalWaveVMinusKey = @"GLInternalWaveVMinusKey";
static NSString *GLInternalWaveWPlusKey = @"GLInternalWaveWPlusKey";
static NSString *GLInternalWaveWMinusKey = @"GLInternalWaveWMinusKey";

@implementation GLInternalWaveInitialization

- (void)encodeWithCoder:(NSCoder *)coder
{
    [coder encodeObject: @(self.maximumModes) forKey:GLInternalWaveMaximumModesKey];
    [coder encodeObject: self.equation forKey:GLInternalWaveEquationKey];
    [coder encodeObject: self.internalModes forKey:GLInternalWaveInternalModeKey];
    [coder encodeObject: self.fullDimensions forKey:GLInternalWaveFullDimensionsKey];
    [coder encodeObject: self.spectralDimensions forKey:GLInternalWaveSpectralDimensionsKey];
    [coder encodeObject: @(self.latitude) forKey:GLInternalWaveLatitudeKey];
    [coder encodeObject: @(self.f0) forKey:GLInternalWaveF0Key];
    [coder encodeObject: self.rho forKey:GLInternalWaveRhoKey];
    [coder encodeObject: self.N2 forKey:GLInternalWaveN2Key];
    [coder encodeObject: self.eigenfrequencies forKey:GLInternalWaveEigenfrequenciesKey];
    [coder encodeObject: self.eigendepths forKey:GLInternalWaveEigendepthsKey];
    [coder encodeObject: self.rossbyRadius forKey:GLInternalWaveRossbyRadiiKey];
    [coder encodeObject: self.S forKey:GLInternalWaveSKey];
    [coder encodeObject: self.Sprime forKey:GLInternalWaveSprimeKey];
    [coder encodeObject: self.zeta_plus forKey:GLInternalWaveZetaPlusKey];
    [coder encodeObject: self.zeta_minus forKey:GLInternalWaveZetaMinusKey];
    [coder encodeObject: self.rho_plus forKey:GLInternalWaveRhoPlusKey];
    [coder encodeObject: self.rho_minus forKey:GLInternalWaveRhoMinusKey];
    [coder encodeObject: self.u_plus forKey:GLInternalWaveUPlusKey];
    [coder encodeObject: self.u_minus forKey:GLInternalWaveUMinusKey];
    [coder encodeObject: self.v_plus forKey:GLInternalWaveVPlusKey];
    [coder encodeObject: self.v_minus forKey:GLInternalWaveVMinusKey];
    [coder encodeObject: self.w_plus forKey:GLInternalWaveWPlusKey];
    [coder encodeObject: self.w_minus forKey:GLInternalWaveWMinusKey];
}

- (id)initWithCoder:(NSCoder *)decoder
{
    if ((self=[super init])) {
        _maximumModes = [[decoder decodeObjectForKey: GLInternalWaveMaximumModesKey] unsignedIntegerValue];
        _equation = [decoder decodeObjectForKey: GLInternalWaveEquationKey];
        _internalModes = [decoder decodeObjectForKey: GLInternalWaveInternalModeKey];
        _fullDimensions = [decoder decodeObjectForKey: GLInternalWaveFullDimensionsKey];
        _spectralDimensions = [decoder decodeObjectForKey: GLInternalWaveSpectralDimensionsKey];
        _latitude = [[decoder decodeObjectForKey: GLInternalWaveLatitudeKey] doubleValue];
        _f0 = [[decoder decodeObjectForKey: GLInternalWaveF0Key] doubleValue];
        _rho = [decoder decodeObjectForKey: GLInternalWaveRhoKey];        
        _N2 = [decoder decodeObjectForKey: GLInternalWaveN2Key];
        _eigenfrequencies = [decoder decodeObjectForKey: GLInternalWaveEigenfrequenciesKey];
        _eigendepths = [decoder decodeObjectForKey: GLInternalWaveEigendepthsKey];
        _rossbyRadius = [decoder decodeObjectForKey: GLInternalWaveRossbyRadiiKey];
        _S = [decoder decodeObjectForKey: GLInternalWaveSKey];
        _Sprime = [decoder decodeObjectForKey: GLInternalWaveSprimeKey];
        _zeta_plus = [decoder decodeObjectForKey: GLInternalWaveZetaPlusKey];
        _zeta_minus = [decoder decodeObjectForKey: GLInternalWaveZetaMinusKey];
        _rho_plus = [decoder decodeObjectForKey: GLInternalWaveRhoPlusKey];
        _rho_minus = [decoder decodeObjectForKey: GLInternalWaveRhoMinusKey];
        _u_plus = [decoder decodeObjectForKey: GLInternalWaveUPlusKey];
        _u_minus = [decoder decodeObjectForKey: GLInternalWaveUMinusKey];
        _v_plus = [decoder decodeObjectForKey: GLInternalWaveVPlusKey];
        _v_minus = [decoder decodeObjectForKey: GLInternalWaveVMinusKey];
        _w_plus = [decoder decodeObjectForKey: GLInternalWaveWPlusKey];
        _w_minus = [decoder decodeObjectForKey: GLInternalWaveWMinusKey];
    }
    return self;
}

- (GLInternalWaveInitialization *) initWithDensityProfile: (GLFunction *) rho fullDimensions: (NSArray *) dimensions latitude: (GLFloat) latitude maxMode: (NSUInteger) maximumModes equation: (GLEquation *) equation
{
	NSUInteger numVerticalDims = 0;
	NSUInteger numHorizontalDims = 0;
	for (GLDimension *dim in dimensions) {
		if ([dim.name isEqualToString: @"x"] || [dim.name isEqualToString: @"y"]) {
			numHorizontalDims++;
		} else if ([dim.name isEqualToString: @"z"]) {
			numVerticalDims++;
		} else {
			[NSException raise: @"BadDimensions" format:@"Dimensions must be name x, y, or z"];
		}
	}
	if (numVerticalDims != 1 || numHorizontalDims == 0 || numHorizontalDims > 2) {
		[NSException raise: @"BadDimensions" format:@"There must be one vertical dimension given, and either 1 or 2 horizontal dimensions."];
	}
	
	GLDimension *zDim = rho.dimensions[0];
	if (maximumModes > zDim.nPoints) {
		maximumModes = zDim.nPoints;
	} else if (maximumModes == 0) {
		if (zDim.nPoints < 8) {
			[NSException raise: @"BadDimensions" format:@"You're going to need 8 or more points in the vertical profile for this to work."];
		} else if (zDim.nPoints > 128) {
			maximumModes = 64;
		} else {
			maximumModes = floor(zDim.nPoints/2);
		}
	}
	
    if ((self=[super init])) {
        self.fullDimensions=dimensions;
        self.equation=equation;
		self.horizontalDimensions = [NSMutableArray array];
		self.rho = rho;
		self.latitude = latitude;
		self.f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
		self.maximumModes = maximumModes;
        
        self.internalModes = [[GLInternalModesSpectral alloc] init];
        [self.internalModes internalWaveModesFromDensityProfile: self.rho withFullDimensions: self.fullDimensions forLatitude: self.latitude maximumModes: maximumModes];
    }
    return self;
}

- (GLInternalWaveInitialization *) initWithDensityProfile: (GLFunction *) rho fullDimensions: (NSArray *) dimensions latitude: (GLFloat) latitude equation: (GLEquation *) equation
{
	return [self initWithDensityProfile: rho fullDimensions: dimensions latitude: latitude maxMode: 0 equation: equation];
}

- (void) computeInternalModes
{
	self.S = self.internalModes.S;
	self.Sprime = self.internalModes.Sprime;
	self.eigenfrequencies = self.internalModes.eigenfrequencies;
	self.eigendepths = self.internalModes.eigendepths;
	self.rossbyRadius = self.internalModes.rossbyRadius;
	self.rho = self.internalModes.rho;
	self.N2 = self.internalModes.N2;
	
	self.rho.name = @"rho_bar";
	self.N2.name = @"N2";
	
	GLFloat minH = [self.eigendepths minNow];
	if (minH <= 0.0) {
		NSLog(@"You have eigendepths with negative values. This is not physical, so you should probably try limiting the number of modes you use.");
	}
	
	GLFloat minTime = 2*M_PI/[self.eigenfrequencies maxNow];
	GLFloat maxTime = 2*M_PI/self.f0;
	NSLog(@"Maximum wave period (based on the Coriolis frequency) is %02d:%02d (HH:MM)", ((int) floor(maxTime/3600))%24, ((int) floor(maxTime/60))%60);
	NSLog(@"Minimum wave period (based on the modes and horizontal resolution) is %02d:%02d (HH:MM)", ((int) floor(minTime/3600))%24, ((int) floor(minTime/60))%60);
	
	self.spectralDimensions = self.eigenfrequencies.dimensions;
	
	for (GLDimension *dim in self.spectralDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			self.kDim = dim;
			[self.horizontalDimensions addObject: dim];
		} else if ( [dim.name isEqualToString: @"l"]) {
			self.lDim = dim;
			[self.horizontalDimensions addObject: dim];
		} else {
			self.modeDim = dim;
		}
	}
}

- (void) showDiagnostics
{
    GLFloat maxOmega = [self.eigenfrequencies maxNow];
    GLFloat dx = self.f0/2;
    NSUInteger N = [self.eigenfrequencies maxNow]/dx;
    
    
    GLDimension *omegaDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints:N+1 domainMin:0 length:N*dx];
    GLFunction *omega = [GLFunction functionOfRealTypeFromDimension: omegaDim withDimensions: @[omegaDim] forEquation: self.equation];
    GLFunction *cardinality = [GLFunction functionOfRealTypeWithDimensions: @[omegaDim] forEquation: self.equation];
    GLFunction *powerU = [GLFunction functionOfRealTypeWithDimensions: @[omegaDim] forEquation: self.equation];
    GLFunction *powerV = [GLFunction functionOfRealTypeWithDimensions: @[omegaDim] forEquation: self.equation];
    [cardinality zero];
    [powerU zero];
    [powerV zero];
    
    GLFunction *up = [self.u_plus abs];
    up = [up multiply: up];
    GLFunction *um = [self.u_minus abs];
    um = [um multiply: um];
    
    GLFunction *vp = [self.v_plus abs];
    vp = [vp multiply: vp];
    GLFunction *vm = [self.v_minus abs];
    vm = [vm multiply: vm];
    
    GLMatrixDescription *md = self.eigenfrequencies.matrixDescription;
    for (NSUInteger mode=0; mode<md.strides[0].nPoints; mode++) {
        for (NSUInteger i=0; i<md.strides[1].nPoints; i++) {
            for (NSUInteger j=0; j<md.strides[2].nPoints; j++) {
                NSUInteger position = mode*md.strides[0].stride + i*md.strides[1].stride + j*md.strides[2].stride;
                NSUInteger bin = self.eigenfrequencies.pointerValue[position]/omegaDim.sampleInterval;
                cardinality.pointerValue[bin] += 1;
                                
                powerU.pointerValue[bin] += up.pointerValue[position];
                powerU.pointerValue[bin] += um.pointerValue[position];
                powerV.pointerValue[bin] += vp.pointerValue[position];
                powerV.pointerValue[bin] += vm.pointerValue[position];
            }
        }
    }
    
    GLFunction *Su = [powerU dividedBy: cardinality];
    GLFunction *Sv = [powerV dividedBy: cardinality];
    
    GLFunction *period = [[omega plus: @(omegaDim.sampleInterval)] scalarDivide: 2*M_PI/60];
    [omega dumpToConsole];
    [cardinality dumpToConsole];
    [Su dumpToConsole];
    [Sv dumpToConsole];
}

- (void) computeInternalModesOld
{
    if (self.maximumModes || self.maxDepth || self.minDepth) {
        // Identify the z dimension
        NSArray *dimensions = self.internalModes.S.toDimensions;
        GLDimension *zDim;
        for (GLDimension *dim in dimensions) {
            if ([dim.name isEqualToString: @"z"]) zDim = dim;
        }
        if (!zDim) {
            [NSException raise:@"BadAssumptions" format:@"Cannot locate a variable name 'z'. This is a requirement."];
        }
        
        NSUInteger maxDepthIndex = 0;
        NSUInteger minDepthIndex = zDim.nPoints-1;;
        
        if (self.maxDepth) {
            for (NSUInteger iPoint = 0; iPoint<zDim.nPoints; iPoint++) {
                if ( [zDim valueAtIndex: iPoint] < self.maxDepth) maxDepthIndex=iPoint;
            }
        }
        
        if (self.minDepth) {
            for (NSInteger iPoint = zDim.nPoints-1; iPoint>=0; iPoint--) {
                if ( [zDim valueAtIndex: iPoint] > self.minDepth) minDepthIndex=iPoint;
            }
        }
        
        NSMutableString *fromIndexString = [NSMutableString stringWithFormat: @""];
        NSMutableString *toIndexString = [NSMutableString stringWithFormat: @""];
        for (GLDimension *dim in dimensions) {
            if (dim==zDim) {
                [fromIndexString appendFormat: @"0:%lu", self.maximumModes-1];
                [toIndexString appendFormat: @"%lu:%lu", maxDepthIndex, minDepthIndex];
            } else {
                [fromIndexString appendFormat: @":"];
                [toIndexString appendFormat: @":"];
            }
            if ([dimensions indexOfObject: dim] < dimensions.count-1) {
                [fromIndexString appendFormat: @","];
                [toIndexString appendFormat: @","];
            }
        }
        self.S = [self.internalModes.S reducedFromDimensions: fromIndexString toDimension: toIndexString];
        self.Sprime = [self.internalModes.Sprime reducedFromDimensions: fromIndexString toDimension: toIndexString];
        self.eigenfrequencies = [self.internalModes.eigenfrequencies variableFromIndexRangeString:fromIndexString];
        self.eigendepths = [self.internalModes.eigendepths variableFromIndexRangeString:fromIndexString];
        self.rossbyRadius = [self.internalModes.rossbyRadius variableFromIndexRangeString:fromIndexString];
        self.rho = [self.internalModes.rho variableFromIndexRangeString:[NSString stringWithFormat: @"%lu:%lu", maxDepthIndex, minDepthIndex]];
        self.N2 = [self.internalModes.N2 variableFromIndexRangeString:[NSString stringWithFormat: @"%lu:%lu", maxDepthIndex, minDepthIndex]];
    } else {
        self.S = self.internalModes.S;
        self.Sprime = self.internalModes.Sprime;
        self.eigenfrequencies = self.internalModes.eigenfrequencies;
        self.eigendepths = self.internalModes.eigendepths;
        self.rossbyRadius = self.internalModes.rossbyRadius;
        self.rho = self.internalModes.rho;
        self.N2 = self.internalModes.N2;
    }
    self.rho.name = @"rho_bar";
    self.N2.name = @"N2";
    
    GLFloat minH = [self.eigendepths minNow];
    if (minH <= 0.0) {
        NSLog(@"You have eigendepths with negative values. This is not physical, so you should probably try limiting the number of modes you use.");
    }
    
    GLFloat minTime = 2*M_PI/[self.eigenfrequencies maxNow];
    GLFloat maxTime = 2*M_PI/self.f0;
    NSLog(@"Maximum wave period (based on the Coriolis frequency) is %02d:%02d (HH:MM)", ((int) floor(maxTime/3600))%24, ((int) floor(maxTime/60))%60);
    NSLog(@"Minimum wave period (based on the modes and horizontal resolution) is %02d:%02d (HH:MM)", ((int) floor(minTime/3600))%24, ((int) floor(minTime/60))%60);
    
	self.spectralDimensions = self.eigenfrequencies.dimensions;
	
	for (GLDimension *dim in self.spectralDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			self.kDim = dim;
			[self.horizontalDimensions addObject: dim];
		} else if ( [dim.name isEqualToString: @"l"]) {
			self.lDim = dim;
			[self.horizontalDimensions addObject: dim];
		} else {
			self.modeDim = dim;
		}
	}
}

- (GLFloat) createUnitWaveWithSpeed: (GLFloat) U_max verticalMode: (NSUInteger) mode k: (NSUInteger) kUnit l: (NSUInteger) lUnit omegaSign: (GLFloat) sign
{
	[self computeInternalModes];
	
    if (mode == 0) {
        [NSException raise: @"UnsupportedMode" format:@"Barotropic mode is not supported. Set the mode to 1 or greater"];
    }
    // The 1st baroclinic mode is in position 0.
    mode = mode-1;
    
    GLFunction *U_mag = [GLFunction functionOfComplexTypeWithDimensions: self.spectralDimensions forEquation: self.equation];
    [U_mag zero];
	
	// My notation is horrible here. kDim and lDim are actually reversed!
    NSUInteger zDimNPoints = [U_mag.dimensions[0] nPoints];
    NSUInteger kDimNPoints = [U_mag.dimensions[1] nPoints];
    NSUInteger lDimNPoints = [U_mag.dimensions[2] nPoints];

    // We will be setting G, the energy density. So we have to convert U_max into that value.
//    GLFloat k = [self.kDim valueAtIndex: kUnit];
//    GLFloat l = [self.lDim valueAtIndex: lUnit];
	GLFloat omega = self.eigenfrequencies.pointerValue[(mode*kDimNPoints+kUnit)*lDimNPoints+lUnit];
    
    GLFloat ratio = (omega*omega-self.f0*self.f0)/(omega*omega+self.f0*self.f0);
    NSLog(@"PE/KE ratio: %f", ratio);
//
//    GLFloat G = U_max*(k*k+l*l)/(k*omega);

    U_mag = [U_mag setValue: 1.0 atIndices: [NSString stringWithFormat:@"%lu,%lu,%lu", mode, kUnit, lUnit]];
	
    GLFloat *C = U_mag.pointerValue;
    
    // index=i*ny*nz+j*nz+k
    // index=(i*ny+j)*nz+k
    // my notation, (z*kDimNPoints+k)*lDimNPoints+l
    for (NSUInteger z=0; z<zDimNPoints; z++) {
        // Hermitian conjugates
        for (NSUInteger i=1; i<kDimNPoints/2; i++) {
            C[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints] = C[(z*kDimNPoints+i)*lDimNPoints];
            C[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)] = C[(z*kDimNPoints+(kDimNPoints-i))*lDimNPoints+(lDimNPoints-1)];
        }
        
        // For the four self-conjugate components, that means that there can be no imaginary component
        // But that their real components should be doubled, to make all else equal.
        C[z*kDimNPoints*lDimNPoints+0] *= 2;
        C[z*kDimNPoints*lDimNPoints+(lDimNPoints-1)] *= 2;
        C[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+0] *= 2;
        C[(z*kDimNPoints+(kDimNPoints/2))*lDimNPoints+(lDimNPoints-1)] *= 2;
    }
    
	//GLScalar *i = [GLScalar scalarWithValue: -I forEquation: self.equation];
	GLFunction *G_plus = [[[U_mag duplicate] makeHermitian] negate];
    GLFunction *G_minus = [[U_mag duplicate] makeHermitian];
    
    if (sign<0) {
        [G_plus zero];
    } else if (sign > 0) {
        [G_minus zero];
    }
    
    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
    
    // Rather than figure out how to properly normalize, we just cheat.
	GLFunction *u = [[self.Sprime transform: [self.u_plus plus: self.u_minus]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
	GLFunction *v = [[self.Sprime transform: [self.v_plus plus: self.v_minus]] transformToBasis: @[@(kGLDeltaBasis), @(kGLDeltaBasis), @(kGLDeltaBasis)]];
	GLFunction *speed = [[[u times: u] plus: [v times: v]] sqrt];
    GLFloat Speed = [speed maxNow];
    
    G_plus = [G_plus times: @(U_max/(Speed))];
    G_minus = [G_minus times: @(U_max/(Speed))];
    
    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
	
	return omega;
}

- (void) createGarrettMunkSpectrumWithEnergy: (GLFloat) energyLevel
{
    [self computeInternalModes];
    
    GLFloat j_star = 3;
    
    GLFloat L_GM = 1.3E3;       // thermocline exponential scale, [m]
    GLFloat invT_GM = 5.2E-3;   // reference buoyancy frequency, [radians/s]
    GLFloat E_GM = 6.3E-5;      // non-dimensional energy parameter, [unitless]
    GLFloat E = (L_GM*L_GM*L_GM*invT_GM*invT_GM*E_GM)*energyLevel; // [ m^3/s^2 ]
    
    // Compute the proper normalization with lots of modes
    // The mode dimension, j, starts at zero, but we want it to start at 1... so we add 1!
    GLDimension *jDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1024 domainMin: 1.0 length:1024.0];
    GLFunction *j1D = [GLFunction functionOfRealTypeFromDimension: jDim withDimensions: @[jDim] forEquation: self.equation];
    GLFunction *H1D = [[[j1D plus: @(j_star)] pow: 5./2.] scalarDivide: 1]; // 3*pow(j_star,3/2)/2
    GLFunction *H_norm = [H1D sum: 0];
    
    j1D = [[GLFunction functionOfRealTypeFromDimension: self.modeDim withDimensions: @[self.modeDim] forEquation: self.equation] plus: @(1)];
    H1D = [[[[j1D plus: @(j_star)] pow: 5./2.] multiply: H_norm] scalarDivide: 1]; // 3*pow(j_star,3/2)/2
    GLScalar *H1D_sum = [H1D sum: 0];
    NSLog(@"The model uses %ld vertical modes, which accounts for %.2f%% of the variance.\n",self.modeDim.nPoints,100.0*H1D_sum.pointerValue[0]);
    
    GLFloat H_norm_scalar = H_norm.pointerValue[0];
    GLFloat f0 = self.f0;
    GLFloat maxOmega = [self.eigenfrequencies maxNow];
    GLFloat B_norm = 1./atan( sqrt(maxOmega*maxOmega/self.f0/self.f0 - 1)); // GM79 assumes this is 2/pi
    GLFloat (^GM2D_omega_function)(GLFloat,GLFloat,GLFloat) = ^(GLFloat omega0, GLFloat omega1, GLFloat j){
        if (omega0 < f0 && omega1 <= f0) {
            return 0.0;
        } else if (omega0 < f0 && omega1 > f0) {
            omega0 = f0; // set the lower limit to the asymptote
        }
        GLFloat B0 = B_norm*atan(f0/sqrt(omega0*omega0 - f0*f0));
        GLFloat B1 = B_norm*atan(f0/sqrt(omega1*omega1 - f0*f0));
        GLFloat H = pow(j+j_star, -5./2.)/H_norm_scalar;
        return -E*H*(B1-B0);
    };
    
    GLFloat maxEnergy = 0.0;
    for (NSUInteger iMode=0; iMode < j1D.nDataPoints; iMode++) {
        maxEnergy += GM2D_omega_function(f0,maxOmega,iMode+1);
    }
    NSLog(@"The GM function sums to a maximum of %f.",maxEnergy/E);
    
    GLFloat maxModeSpeed = 9.81 * [self.eigendepths maxNow];
    GLFloat dk = self.kDim.sampleInterval;
    GLFloat domega = sqrt( maxModeSpeed*dk*dk + f0*f0 );
    
    GLDimension *omegaDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints:1+(maxOmega-f0)/domega domainMin:f0 length:maxOmega];
    GLFunction *omega = [GLFunction functionOfRealTypeFromDimension: omegaDim withDimensions: @[omegaDim] forEquation: self.equation];
    GLFunction *cardinality = [GLFunction functionOfRealTypeWithDimensions: @[omegaDim] forEquation: self.equation];
    
    // GM3D stores the *expected* energy of each wave for a given (j,k,l)
    GLFunction *GM3D = [GLFunction functionOfRealTypeWithDimensions: self.spectralDimensions forEquation: self.equation];
    
    GLMatrixDescription *md = self.eigenfrequencies.matrixDescription;
    maxEnergy = 0.0;
    for (NSUInteger mode=0; mode<md.strides[0].nPoints; mode++) {
        
        // Loop over our omega bins
        for (NSUInteger w=0; w<omegaDim.nPoints-1; w++ ) {
            GLFloat omega_lower = omega.pointerValue[w];
            GLFloat omega_upper = omega.pointerValue[w+1];
            
            // Find out how many frequencies are in this band
            NSUInteger totalFrequenciesInRange = 0;
            for (NSUInteger i=0; i<md.strides[1].nPoints; i++) {
                for (NSUInteger j=0; j<md.strides[2].nPoints; j++) {
                    NSUInteger position = mode*md.strides[0].stride + i*md.strides[1].stride + j*md.strides[2].stride;
                    if (self.eigenfrequencies.pointerValue[position] >= omega_lower && self.eigenfrequencies.pointerValue[position] < omega_upper) {
                        totalFrequenciesInRange += j == 0 ? 1 : 2;
                    }
                }
            }
            maxEnergy += GM2D_omega_function(omega_lower,omega_upper,mode+1);
            GLFloat energyPerFrequency = GM2D_omega_function(omega_lower,omega_upper,mode+1)/totalFrequenciesInRange;
            
            // Now put the energy in this band
            for (NSUInteger i=0; i<md.strides[1].nPoints; i++) {
                for (NSUInteger j=0; j<md.strides[2].nPoints; j++) {
                    NSUInteger position = mode*md.strides[0].stride + i*md.strides[1].stride + j*md.strides[2].stride;
                    if (self.eigenfrequencies.pointerValue[position] >= omega_lower && self.eigenfrequencies.pointerValue[position] < omega_upper) {
                        GM3D.pointerValue[position] = energyPerFrequency;
                    }
                }
            }
            
        }
    }
    NSLog(@"The GM function sums to a maximum of %f.",maxEnergy/E);
    NSLog(@"Due to horizontal domain size limitations, %.2f%% of the energy is in the j=1, k=l=0 mode.",100*GM3D.pointerValue[0]/E);

    
    // This will randomize the energy of phase of each wave. The *mean* variance (energy) is 1.0.
    GLFunction *G_plus = [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation];
    GLFunction *G_minus = [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation];
    
    NSUInteger numReps = 1;
    for (NSUInteger i=0; i<numReps-1; i++) {
        G_plus = [G_plus plus: [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation]];
        G_minus = [G_minus plus: [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation]];
    }
    G_plus = [G_plus times: @(sqrt(1/((GLFloat)numReps)))];
    G_minus = [G_minus times: @(sqrt(1/((GLFloat)numReps)))];
    
    // This forces amplitudes to be uniform, but phases to remain randomized. Good for diagnostics.
    G_plus = [G_plus dividedBy: [G_plus abs]];
    G_minus = [G_minus dividedBy: [G_minus abs]];
    
    
    
    GLFunction *GM_sum = [[[[GM3D times: @(1/E)] sum: 2] sum: 1] sum: 0];
    GLFunction *Conjugacy = [GM3D variableFromIndexRangeString: @":,:,1:end"];
    GLFunction *GM_Conjugacy_sum = [[[[Conjugacy times: @(1/E)] sum: 2] sum: 1] sum: 0];
    
    NSLog(@"The GM coefficients should sum to 1. In practice, they sum to: %.2f%%.",100.*(GM_sum.pointerValue[0]+GM_Conjugacy_sum.pointerValue[0]));
    
    // We cut the value in half, because we will want the expectation of the the positive and negative sides to add up to this value.
    GLFunction *G = [[GM3D times: @(0.5)] sqrt] ;
    
    [self zeroOutUnphysicalComponentsInAmplitudeFunction: G];
    
    // We seed the random number generator, then call and resolve the randomized function so that their order isn't flipped.
    srand(4);
    
    
    
    GLFunction *G_sum1 = [[[[[G_plus abs] pow: 2.0] sum: 2] sum: 1] sum: 0];
    GLFunction *G_sum2 = [[[[[G_minus abs] pow: 2.0] mean: 2] mean: 1] mean: 0];
    NSLog(@"The random coefficients sum to: %f, %f",G_sum1.pointerValue[0],G_sum2.pointerValue[0]);
    
    G_plus = [G_plus multiply: G];
    G_minus = [G_minus multiply: G];
    
    G_sum1 = [[[[[G_plus abs] pow: 2.0] sum: 2] sum: 1] sum: 0]; GLFunction *G_sum1_conjugacy = [[[[[[G_plus abs] pow: 2.0] variableFromIndexRangeString: @"1:end,1:end,:"]sum: 2] sum: 1] sum: 0];
    G_sum2 = [[[[[G_minus abs] pow: 2.0] sum: 2] sum: 1] sum: 0]; GLFunction *G_sum2_conjugacy = [[[[[[G_minus abs] pow: 2.0] variableFromIndexRangeString: @"1:end,1:end,:"]sum: 2] sum: 1] sum: 0];
    NSLog(@"The random coefficients sum to: %f, %f",(G_sum1.pointerValue[0]+G_sum1_conjugacy.pointerValue[0])/E,(G_sum2.pointerValue[0]+G_sum2_conjugacy.pointerValue[0])/E);
    
    GLFunction * GM_random_sum = [[[[[[[G_plus abs] multiply: [G_plus abs]] plus: [[G_minus abs] multiply: [G_minus abs]]] times: @(1/E)] sum: 2] sum: 1] sum: 0];
    NSLog(@"The GM random coefficients should sum to a value near 1. In practice, they sum to: %f",GM_random_sum.pointerValue[0]);
    
    //G_minus = [G_minus setValue: 0.0 atIndices: @":,0,0"]; // Inertial motions go only one direction!
    
    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
    
    
    printf("\n\n[");
    domega = 2*M_PI/(0.5*86400);
    NSUInteger t=0;
    GLFloat totalEnergyInFrequencyDomain = 0.0;
    for (GLFloat omega=0.0; omega <= maxOmega; omega += domega) {
        GLFloat totalEnergyInWaveband = 0.0;
        for (NSUInteger iMode=0; iMode < j1D.nDataPoints; iMode++) {
            for (NSUInteger i=0; i<md.strides[1].nPoints; i++) {
                for (NSUInteger j=0; j<md.strides[2].nPoints; j++) {
                    NSUInteger index = iMode*md.strides[0].stride + i*md.strides[1].stride + j*md.strides[2].stride;
                    if (self.eigenfrequencies.pointerValue[index] >= omega && self.eigenfrequencies.pointerValue[index] < omega+domega) {
                        totalEnergyInWaveband += (j!=0 ? 2.0 : 1.0)*GM3D.pointerValue[index];
                        t++;
                    }
                }
            }
        }
        totalEnergyInFrequencyDomain += totalEnergyInWaveband;
        printf("%f,%g;\n",omega,totalEnergyInWaveband);
    }
    printf("];\n\n");
    
    NSLog(@"Visited %lu of %lu points. Summed the energy to %.2f%% of GM.",t,GM3D.nDataElements,100*totalEnergyInFrequencyDomain/E);
    
}

//- (void) createGarrettMunkSpectrumWithEnergy: (GLFloat) energyLevel
//{
//	[self computeInternalModes];
//	
//    GLFloat j_star = 3;
//    
//    GLFloat L_GM = 1.3E3;       // thermocline exponential scale, [m]
//    GLFloat invT_GM = 5.2E-3;   // reference buoyancy frequency, [radians/s]
//    GLFloat E_GM = 6.3E-5;      // non-dimensional energy parameter, [unitless]
//    GLFloat E = (L_GM*L_GM*L_GM*invT_GM*invT_GM*E_GM)*energyLevel; // [ m^3/s^2 ]
//	
//    // Compute the proper normalization with lots of modes
//    // The mode dimension, j, starts at zero, but we want it to start at 1... so we add 1!
//    GLDimension *jDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 1024 domainMin: 1.0 length:1024.0];
//    GLFunction *j1D = [GLFunction functionOfRealTypeFromDimension: jDim withDimensions: @[jDim] forEquation: self.equation];
//    GLFunction *H1D = [[[j1D plus: @(j_star)] pow: 5./2.] scalarDivide: 1]; // 3*pow(j_star,3/2)/2
//    GLFunction *H_norm = [H1D sum: 0];
//    
//    j1D = [[GLFunction functionOfRealTypeFromDimension: self.modeDim withDimensions: @[self.modeDim] forEquation: self.equation] plus: @(1)];
//    H1D = [[[[j1D plus: @(j_star)] pow: 5./2.] multiply: H_norm] scalarDivide: 1]; // 3*pow(j_star,3/2)/2
//    GLScalar *H1D_sum = [H1D sum: 0];
//    NSLog(@"The model uses %ld vertical modes, which accounts for %.2f%% of the variance.\n",self.modeDim.nPoints,100.0*H1D_sum.pointerValue[0]);
//    
//    
//    // Compute the total energy between two wavenumbers (k0, k1) for a given mode j.
//    GLScalar *intN = [[self.N2 sqrt] integrateToLimits];
//    GLScalar *intInvN = [[[self.N2 sqrt] scalarDivide: 1.0] integrateToLimits];
//    GLScalar *e2 = [intN dividedBy: intInvN];
//    GLScalar *s2 = [intN times: intInvN];
//    GLFloat eta2 = e2.pointerValue[0];
//    GLFloat maxOmega = [self.eigenfrequencies maxNow];
//    NSLog(@"Your WKB buoyancy frequency is %.2f cycles per hour, actual is %.2f.", sqrt(eta2)*3600/(2*pi), maxOmega*3600/(2*pi));
//    
//    GLFloat m_norm = M_PI*M_PI/s2.pointerValue[0];
//    GLFloat H_norm_scalar = H_norm.pointerValue[0];
//    GLFloat f0 = self.f0;
//    GLFloat B_norm = 1./atan( sqrt(eta2/self.f0/self.f0 - 1)); // GM79 assumes this is 2/pi
//    GLFloat (^GM2D_function)(GLFloat,GLFloat,GLFloat) = ^(GLFloat k0, GLFloat k1, GLFloat j){
//        GLFloat mj2 = j*j*m_norm;
//        GLFloat B0 = B_norm*atan((k0/f0)*sqrt((eta2 - f0*f0)/(k0*k0+mj2)));
//        GLFloat B1 = B_norm*atan((k1/f0)*sqrt((eta2 - f0*f0)/(k1*k1+mj2)));
//        GLFloat H = pow(j+j_star, -5./2.)/H_norm_scalar;
//        return E*H*(B1-B0);
//    };
//    
//    GLFloat (^GM2D_omega_function)(GLFloat,GLFloat,GLFloat) = ^(GLFloat omega0, GLFloat omega1, GLFloat j){
//        if (omega0 < f0 && omega1 <= f0) {
//            return 0.0;
//        } else if (omega0 < f0 && omega1 > f0) {
//            omega0 = f0; // set the lower limit to the asymptote
//        }
//        GLFloat mj2 = j*j*m_norm;
//        GLFloat B0 = B_norm*atan(f0/sqrt(omega0*omega0 - f0*f0));
//        GLFloat B1 = B_norm*atan(f0/sqrt(omega1*omega1 - f0*f0));
//        GLFloat H = pow(j+j_star, -5./2.)/H_norm_scalar;
//        return E*H*(B1-B0);
//    };
//    
//    GLFloat maxEnergy = 0.0;
//    for (NSUInteger iMode=0; iMode < j1D.nDataPoints; iMode++) {
//        maxEnergy += GM2D_function(0,10000*2*M_PI*self.kDim.domainLength,iMode+1);
//    }
//    NSLog(@"The GM function sums to a maximum of %f. If this is less than 1, the reason for this is that B(omega) is not normalized correctly. Integral needs to go to N, not infinity.",maxEnergy/E);
//    
//	GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: self.kDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
//	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: self.lDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
//    GLFunction *Kh = [[[k multiply: k] plus: [l multiply: l]] sqrt];
//	
//    // GM3D stores the *expected* energy of each wave for a given (j,k,l)
//    GLFunction *GM3D = [GLFunction functionOfRealTypeWithDimensions: self.spectralDimensions forEquation: self.equation];
//    
//    // This will randomize the energy of phase of each wave. The *mean* variance (energy) is 1.0.
//    GLFunction *G_plus = [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation];
//    GLFunction *G_minus = [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation];
//    
//    NSUInteger numReps = 1;
//    for (NSUInteger i=0; i<numReps-1; i++) {
//        G_plus = [G_plus plus: [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation]];
//        G_minus = [G_minus plus: [GLFunction functionWithNormallyDistributedValueWithDimensions: self.spectralDimensions forEquation: self.equation]];
//    }
//    G_plus = [G_plus times: @(sqrt(1/((GLFloat)numReps)))];
//    G_minus = [G_minus times: @(sqrt(1/((GLFloat)numReps)))];
//    
//    // This forces amplitudes to be uniform, but phases to remain randomized. Good for diagnostics.
//    G_plus = [G_plus dividedBy: [G_plus abs]];
//    G_minus = [G_minus dividedBy: [G_minus abs]];
//    
//    // 3D real FFTs are hard b/c approximately half the energy is assumed from Hermitian conjugacy.
//    // using (i,j,k) ordering for dims (nx,ny,nz), note that the points with assumed (not explicit) conjugacy are:
//    // page/row/column -- for each page, columns > 0
//    NSUInteger jDimNPoints = [GM3D.dimensions[0] nPoints];
//    NSUInteger kDimNPoints = [GM3D.dimensions[1] nPoints];
//    NSUInteger lDimNPoints = [GM3D.dimensions[2] nPoints];
//    GLFloat missingEnergy = 0.0; GLFloat max_m = 0.0;
//    GLFloat notMissingEnergy = 0.0;
//    for (NSUInteger iMode=0; iMode < j1D.nDataPoints; iMode++) {
//        // Set k=l=0 wavenumber for this mode
//        GM3D.pointerValue[iMode*kDimNPoints*lDimNPoints] = GM2D_function(0,2*M_PI*self.kDim.sampleInterval/2,iMode+1);
//        notMissingEnergy += GM2D_function(0,2*M_PI*self.kDim.sampleInterval/2,iMode+1);
//        
//        // Walk through the other wavenumbers in a HORRIBLY inefficient fashion and distribute the energy;
//        for (NSUInteger m=1; m<self.kDim.nPoints/2; m++ ) {
//            GLFloat m_lower = 2*M_PI*([self.kDim valueAtIndex: m] - self.kDim.sampleInterval/2);
//            GLFloat m_upper = 2*M_PI*([self.kDim valueAtIndex: m] + self.kDim.sampleInterval/2);
//            
//            NSUInteger totalWavenumbersInRange = 0;
//            NSUInteger totalRandomRepsInRange = 0;
//            for (NSUInteger i=0; i<kDimNPoints; i++) {
//                for (NSUInteger j=0; j<lDimNPoints; j++) {
//                    NSUInteger index = (iMode*kDimNPoints+i)*lDimNPoints+j; // (i*ny+j)*nz+k
//                    if (Kh.pointerValue[index] >= m_lower && Kh.pointerValue[index] < m_upper) {
//                        totalWavenumbersInRange += j == 0 ? 1 : 2; // This accounts for the hermitian conjugacy energy
//                        totalRandomRepsInRange += 1;
//                    }
//                }
//            }
//            
//            GLFloat energyPerWavenumber = GM2D_function(m_lower,m_upper,iMode+1)/totalWavenumbersInRange;
//            notMissingEnergy += GM2D_function(m_lower,m_upper,iMode+1);
//            
//            for (NSUInteger i=0; i<kDimNPoints; i++) {
//                for (NSUInteger j=0; j<lDimNPoints; j++) {
//                    NSUInteger index = (iMode*kDimNPoints+i)*lDimNPoints+j; // (i*ny+j)*nz+k
//                    if (Kh.pointerValue[index] >= m_lower && Kh.pointerValue[index] < m_upper) {
//                        GM3D.pointerValue[index] = energyPerWavenumber;
//                        
////                        NSUInteger numRepsRequired = floor(100.0*energyPerWavenumber/E);
////                        if (numRepsRequired > 0)
////                        {   // If the energy in this wavenumber band exceeds 1%, reduce the variance of the random number, so it doesn't too heavily effect total energy.
////                            for (NSUInteger rep=0; rep<numRepsRequired; rep++) {
////                                
////                            }
////                        }
//                    }
//                }
//            }
//            
//            max_m = m_upper;
//        }
//        
//        missingEnergy += GM2D_function(max_m,1000*max_m,iMode+1);
//    }
//    
//    // Note to self --- things do add up correctly when we go to large domain sizes.
//    NSLog(@"Due to horizontal grid resolution, %.2f%% of the energy is missing (%.2f%% is accounted for).",100*missingEnergy/E,100*notMissingEnergy/E);
//    NSLog(@"Due to horizontal domain size limitations, %.2f%% of the energy is in the j=k=l=0 mode.",100*GM3D.pointerValue[0]/E);
//    
//    GLFunction *GM_sum = [[[[GM3D times: @(1/E)] sum: 2] sum: 1] sum: 0];
//    GLFunction *Conjugacy = [GM3D variableFromIndexRangeString: @":,:,1:end"];
//    GLFunction *GM_Conjugacy_sum = [[[[Conjugacy times: @(1/E)] sum: 2] sum: 1] sum: 0];
//    
//    NSLog(@"The GM coefficients should sum to 1. In practice, they sum to: %.2f%%.",100.*(GM_sum.pointerValue[0]+GM_Conjugacy_sum.pointerValue[0]));
//    
//	// We cut the value in half, because we will want the expectation of the the positive and negative sides to add up to this value.
//	GLFunction *G = [[GM3D times: @(0.5)] sqrt] ;
//	
//    [self zeroOutUnphysicalComponentsInAmplitudeFunction: G];
//    
//    // We seed the random number generator, then call and resolve the randomized function so that their order isn't flipped.
//    srand(4);
//    
//    
//    
//    GLFunction *G_sum1 = [[[[[G_plus abs] pow: 2.0] sum: 2] sum: 1] sum: 0];
//    GLFunction *G_sum2 = [[[[[G_minus abs] pow: 2.0] mean: 2] mean: 1] mean: 0];
//    NSLog(@"The random coefficients sum to: %f, %f",G_sum1.pointerValue[0],G_sum2.pointerValue[0]);
//    
//	G_plus = [G_plus multiply: G];
//	G_minus = [G_minus multiply: G];
//    
//    G_sum1 = [[[[[G_plus abs] pow: 2.0] sum: 2] sum: 1] sum: 0]; GLFunction *G_sum1_conjugacy = [[[[[[G_plus abs] pow: 2.0] variableFromIndexRangeString: @"1:end,1:end,:"]sum: 2] sum: 1] sum: 0];
//    G_sum2 = [[[[[G_minus abs] pow: 2.0] sum: 2] sum: 1] sum: 0]; GLFunction *G_sum2_conjugacy = [[[[[[G_minus abs] pow: 2.0] variableFromIndexRangeString: @"1:end,1:end,:"]sum: 2] sum: 1] sum: 0];
//    NSLog(@"The random coefficients sum to: %f, %f",(G_sum1.pointerValue[0]+G_sum1_conjugacy.pointerValue[0])/E,(G_sum2.pointerValue[0]+G_sum2_conjugacy.pointerValue[0])/E);
//    
//    GLFunction * GM_random_sum = [[[[[[[G_plus abs] multiply: [G_plus abs]] plus: [[G_minus abs] multiply: [G_minus abs]]] times: @(1/E)] sum: 2] sum: 1] sum: 0];
//    NSLog(@"The GM random coefficients should sum to a value near 1. In practice, they sum to: %f",GM_random_sum.pointerValue[0]);
//          
//    //G_minus = [G_minus setValue: 0.0 atIndices: @":,0,0"]; // Inertial motions go only one direction!
//	
//    [self generateWavePhasesFromPositive: G_plus negative: G_minus];
//    
//    
//    printf("\n\n[");
//    GLFloat domega = 2*M_PI/(0.5*86400);
//    NSUInteger t=0;
//    GLFloat totalEnergyInFrequencyDomain = 0.0;
//    for (GLFloat omega=0.0; omega <= maxOmega; omega += domega) {
//        GLFloat totalEnergyInWaveband = 0.0;
//        for (NSUInteger iMode=0; iMode < j1D.nDataPoints; iMode++) {
//            for (NSUInteger i=0; i<kDimNPoints; i++) {
//                for (NSUInteger j=0; j<lDimNPoints; j++) {
//                    NSUInteger index = (iMode*kDimNPoints+i)*lDimNPoints+j; // (i*ny+j)*nz+k
//                    if (self.eigenfrequencies.pointerValue[index] >= omega && self.eigenfrequencies.pointerValue[index] < omega+domega) {
//                        totalEnergyInWaveband += (j!=0 ? 2.0 : 1.0)*GM3D.pointerValue[index];
//                        t++;
//                    }
//                }
//            }
//        }
//        totalEnergyInFrequencyDomain += totalEnergyInWaveband;
//        printf("%f,%g;\n",omega,totalEnergyInWaveband);
//    }
//    printf("];\n\n");
//    
//    NSLog(@"Visited %lu of %lu points. Summed the energy to %.2f%% of GM.",t,GM3D.nDataElements,100*totalEnergyInFrequencyDomain/E);
//    
//}

- (void) zeroOutUnphysicalComponentsInAmplitudeFunction: (GLFunction *) G
{
    NSUInteger zeroed = 0;
    NSUInteger notZeroed = 0;
    GLFloat Nmax = sqrt(fabs([self.N2 maxNow]));
    for (NSUInteger i=0; i<self.eigenfrequencies.nDataPoints; i++) {
        if (self.eigenfrequencies.pointerValue[i] > Nmax) {
            G.pointerValue[i] = 0;
            zeroed++;
        } else {
            notZeroed++;
        }
    }
    NSLog(@"Zeroed the amplitude %lu frequencies, left %lu untouched because these components have a frequency greater than the buoyancy frequency.", zeroed, notZeroed);
}

// G_plus and G_minus must have Hermitian symmetry, otherwise this won't work.
- (void) generateWavePhasesFromPositive: (GLFunction *) G_plus negative: (GLFunction *) G_minus
{
    GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: self.kDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: self.lDim withDimensions: self.spectralDimensions forEquation:self.equation] scalarMultiply: 2*M_PI];
    GLFunction *K_H = [[[k multiply: k] plus: [l multiply: l]] sqrt];
    GLFunction *sqrtH = [self.eigendepths sqrt];
    
	GLScalar *rho0 = [self.rho min];
    
	self.zeta_plus = [G_plus multiply: [[[K_H multiply: sqrtH] dividedBy: self.eigenfrequencies] makeHermitian]];
    self.zeta_minus = [G_minus multiply: [[[K_H multiply: sqrtH] dividedBy: self.eigenfrequencies] makeHermitian]];
	
	self.rho_plus = [[self.zeta_plus times: rho0] times: @(1/g)];
    self.rho_minus = [[self.zeta_minus times: rho0] times: @(1/g)];
	
	self.w_plus = [G_plus multiply: [[[K_H multiply: sqrtH] swapComplex] makeHermitian]];
    self.w_minus = [G_minus multiply: [[[[K_H multiply: sqrtH] swapComplex] negate] makeHermitian]];
    
    GLFunction *alpha = [l atan2: k];
    GLFunction *cosAlpha = [alpha cos];
    GLFunction *sinAlpha = [alpha sin];
    GLFunction *denominator = [self.eigenfrequencies multiply: sqrtH];
    
	self.u_plus = [G_plus multiply: [[[[[cosAlpha multiply: self.eigenfrequencies] minus: [[sinAlpha times: @(self.f0)] swapComplex]] dividedBy: denominator] negate] makeHermitian]];
    self.u_minus = [G_minus multiply: [[[[cosAlpha multiply: self.eigenfrequencies] plus: [[sinAlpha times: @(self.f0)] swapComplex]] dividedBy: denominator] makeHermitian]];
	self.u_minus = [self.u_minus setValue: 0 atIndices: @":,0,0"]; // Inertial motions go only one direction! This is a special case of the solution.
    
	self.v_plus = [G_plus multiply: [[[[[sinAlpha multiply: self.eigenfrequencies] plus: [[cosAlpha times: @(self.f0)] swapComplex]] dividedBy: denominator] negate] makeHermitian]];
    self.v_minus = [G_minus multiply: [[[[sinAlpha multiply: self.eigenfrequencies] minus: [[cosAlpha times: @(self.f0)] swapComplex]] dividedBy: denominator] makeHermitian]];
#warning Need to cancel out v_minus!!!
    
    self.zeta_plus.name = @"zeta_plus";
    self.zeta_minus.name = @"zeta_minus";
    self.rho_plus.name = @"rho_plus";
    self.rho_minus.name = @"rho_minus";
    self.w_plus.name = @"w_plus";
    self.w_minus.name = @"w_minus";
    self.u_plus.name = @"u_plus";
    self.u_minus.name = @"u_minus";
    self.v_plus.name = @"v_plus";
    self.v_minus.name = @"v_minus";
}


@end
