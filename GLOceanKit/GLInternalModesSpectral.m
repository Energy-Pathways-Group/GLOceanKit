//
//  GLInternalModesSpectral.m
//  GLOceanKit
//
//  Created by Jeffrey J. Early on 3/5/15.
//  Copyright (c) 2015 Jeffrey J. Early. All rights reserved.
//

#import "GLInternalModesSpectral.h"

#define g 9.81

@interface GLInternalModesSpectral ()
- (void) createStratificationProfileFromDensity: (GLFunction *) rho atLatitude: (GLFloat) latitude;
- (void) normalizeDepthBasedEigenvalues: (GLFunction *) lambda eigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm;
- (void) normalizeFrequencyBasedEigenvalues: (GLFunction *) lambda eigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm;
- (void) normalizeEigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm;
@property(strong) GLEquation *equation;

/// The dimension the user used for rho
@property(strong) GLDimension *zInDim;

/// A ChebyshevEndpointGrid dimension similar to the users input dimension.
@property(strong) GLDimension *zDim;

/// Transform of the the above zDim
@property(strong) GLDimension *chebDim;

/// A copy of rho projected (interpolated) onto zDim
@property(strong) GLFunction *rho_cheb_grid;

/// N2 on zDim
@property(strong) GLFunction *N2_cheb_grid;

// These are all on zDim
@property(strong) GLLinearTransform *diffZ;
@property(strong) GLLinearTransform *T;
@property(strong) GLLinearTransform *Tz;
@property(strong) GLLinearTransform *Tzz;

// Chebyshev polynomials on the zInDim.
@property(strong) GLLinearTransform *T_zIn;

// The 'expanded' dimensions.
@property(strong) NSArray *fromDimensions;
@property(strong) NSArray *toDimensions;
@property(strong) NSArray *toFinalDimensions;


@property(strong) GLOperationOptimizer *optimizer;
@property(copy) variableOperation operation;
@property(strong) NSMutableArray *dataBuffers;
@property(strong) GLScalar *k2;

@end

@implementation GLInternalModesSpectral

static NSString *GLInternalModeF0Key = @"GLInternalModeF0Key";
static NSString *GLInternalModeRhoKey = @"GLInternalModeRhoKey";
static NSString *GLInternalModeN2Key = @"GLInternalModeN2Key";
static NSString *GLInternalModeEigenfrequenciesKey = @"GLInternalModeEigenfrequenciesKey";
static NSString *GLInternalModeEigendepthsKey = @"GLInternalModeEigendepthsKey";
static NSString *GLInternalModeRossbyRadiiKey = @"GLInternalModeRossbyRadiiKey";
static NSString *GLInternalModeSKey = @"GLInternalModeSKey";
static NSString *GLInternalModeSprimeKey = @"GLInternalModeSprimeKey";
static NSString *GLInternalModeKDimKey = @"GLInternalModeKDimKey";
static NSString *GLInternalModeLDimKey = @"GLInternalModeLDimKey";

- (void)encodeWithCoder:(NSCoder *)coder
{
	[coder encodeObject: @(self.f0) forKey:GLInternalModeF0Key];
	[coder encodeObject: self.rho forKey:GLInternalModeRhoKey];
	[coder encodeObject: self.N2 forKey:GLInternalModeN2Key];
	[coder encodeObject: self.eigenfrequencies forKey:GLInternalModeEigenfrequenciesKey];
	[coder encodeObject: self.eigendepths forKey:GLInternalModeEigendepthsKey];
	[coder encodeObject: self.rossbyRadius forKey:GLInternalModeRossbyRadiiKey];
	[coder encodeObject: self.S forKey:GLInternalModeSKey];
	[coder encodeObject: self.Sprime forKey:GLInternalModeSprimeKey];
	[coder encodeObject: self.k forKey:GLInternalModeKDimKey];
	[coder encodeObject: self.l forKey:GLInternalModeLDimKey];
}

- (id)initWithCoder:(NSCoder *)decoder
{
	if ((self=[super init])) {
		_f0 = [[decoder decodeObjectForKey: GLInternalModeF0Key] doubleValue];
		_rho = [decoder decodeObjectForKey: GLInternalModeRhoKey];
		_N2 = [decoder decodeObjectForKey: GLInternalModeN2Key];
		_eigenfrequencies = [decoder decodeObjectForKey: GLInternalModeEigenfrequenciesKey];
		_eigendepths = [decoder decodeObjectForKey: GLInternalModeEigendepthsKey];
		_rossbyRadius = [decoder decodeObjectForKey: GLInternalModeRossbyRadiiKey];
		_S = [decoder decodeObjectForKey: GLInternalModeSKey];
		_Sprime = [decoder decodeObjectForKey: GLInternalModeSprimeKey];
		_k = [decoder decodeObjectForKey: GLInternalModeKDimKey];
		_l = [decoder decodeObjectForKey: GLInternalModeLDimKey];
	}
	return self;
}

- (void) createStratificationProfileFromDensity: (GLFunction *) rho atLatitude: (GLFloat) latitude
{
	if (rho.dimensions.count != 1) {
		[NSException raise:@"InvalidDimensions" format:@"Only one dimension allowed, at this point"];
	}
	
	GLScalar *rho0 = [rho min];
	self.f0 = 2*(7.2921e-5)*sin(latitude*M_PI/180);
	self.equation = rho.equation;
	
	GLDimension *zInitial = rho.dimensions[0];
	if (zInitial.gridType == kGLChebyshevEndpointGrid) {
        self.zInDim = zInitial;
		self.zDim = zInitial;
		self.rho = rho;
        self.rho_cheb_grid = rho;
	} else {
        self.zInDim = zInitial;
		self.zDim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: zInitial.nPoints domainMin:zInitial.domainMin length:zInitial.domainLength];
		GLFunction *z = [GLFunction functionOfRealTypeFromDimension: self.zDim withDimensions: @[self.zDim] forEquation: self.equation];
        self.rho = rho;
		self.rho_cheb_grid = [rho interpolateAtPoints: @[z]];
	}
    
    GLFunction *rho_cheb = [[self.rho_cheb_grid minus: rho0] transformToBasis: @[@(kGLChebyshevBasis)]];
	GLFunction *buoyancy = [[rho_cheb dividedBy: rho0] times: @(-g)];
	self.chebDim = rho_cheb.dimensions[0];
	
	self.diffZ = [GLLinearTransform differentialOperatorWithDerivatives: 1 fromDimension: self.chebDim forEquation: self.equation];
    GLFunction *buoyancy_z = [self.diffZ transform: buoyancy];
    self.N2_cheb_grid = [buoyancy_z transformToBasis: @[@(kGLDeltaBasis)]];
    
    self.T_zIn = [GLLinearTransform discreteTransformFromDimension: self.chebDim toDimension: self.zInDim forEquation: self.equation];
    
	self.N2 = [self.T_zIn transform: buoyancy_z];
	self.N2.name = @"N2";
    
    [self createChebyshevTransformations];
}

- (void) createChebyshevTransformations
{
    self.T = [GLLinearTransform discreteTransformFromDimension: self.chebDim toBasis: kGLDeltaBasis forEquation:self.equation];
    // Create an operation to compute the Chebyshev differentiation matrices.
    GLMatrixDescription *matrixDescription = self.diffZ.matrixDescription;
    GLFloat scale = 2./self.zDim.domainLength;
    variableOperation chebOp = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        GLFloat *A = (GLFloat *) [operandArray[0] bytes];
        GLFloat *C = (GLFloat *) [resultArray[0] bytes];
        
        NSUInteger j=0;
        for (NSUInteger i=0; i<matrixDescription.strides[0].nRows; i++) {
            C[i*matrixDescription.strides[0].rowStride + j*matrixDescription.strides[0].columnStride] = 0;
        }
        j=1;
        for (NSUInteger i=0; i<matrixDescription.strides[0].nRows; i++) {
            C[i*matrixDescription.strides[0].rowStride + j*matrixDescription.strides[0].columnStride] = 2*A[i*matrixDescription.strides[0].rowStride + (j-1)*matrixDescription.strides[0].columnStride];
        }
        j=2;
        for (NSUInteger i=0; i<matrixDescription.strides[0].nRows; i++) {
            C[i*matrixDescription.strides[0].rowStride + j*matrixDescription.strides[0].columnStride] = 4*A[i*matrixDescription.strides[0].rowStride + (j-1)*matrixDescription.strides[0].columnStride];
        }
        
        for (j=3; j<matrixDescription.strides[0].nColumns; j++) {
            for (NSUInteger i=0; i<matrixDescription.strides[0].nRows; i++) {
                GLFloat coeff = ((GLFloat) j)/(((GLFloat) j)-2.);
                C[i*matrixDescription.strides[0].rowStride + j*matrixDescription.strides[0].columnStride] = coeff*C[i*matrixDescription.strides[0].rowStride + (j-2)*matrixDescription.strides[0].columnStride] + 2*j*A[i*matrixDescription.strides[0].rowStride + (j-1)*matrixDescription.strides[0].columnStride];
            }
        }
        
        vGL_vsmul(C, matrixDescription.strides[0].stride, &scale,C, matrixDescription.strides[0].stride,matrixDescription.strides[0].nPoints);
    };
    
    GLVariableOperation *operation = [[GLVariableOperation alloc] initWithResult: @[[GLVariable variableWithPrototype: self.T]] operand:@[self.T] buffers:@[] operation: chebOp];
    self.Tz = operation.result[0];
    
    GLVariableOperation *op2 = [[GLVariableOperation alloc] initWithResult: @[[GLVariable variableWithPrototype: self.Tz]] operand:@[self.Tz] buffers:@[] operation: chebOp];
    self.Tzz = op2.result[0];
}

- (void) normalizeDepthBasedEigenvalues: (GLFunction *) lambda eigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm
{
	self.eigendepths = [lambda makeRealIfPossible]; self.eigendepths.name = @"eigendepths";
	self.rossbyRadius = [[[self.eigendepths times: @(g/(self.f0*self.f0))] abs] sqrt]; self.rossbyRadius.name = @"rossbyRadii";
	[self normalizeEigenvectors:S withNorm: norm];
}

- (void) normalizeFrequencyBasedEigenvalues: (GLFunction *) lambda eigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm
{
	lambda = [lambda makeRealIfPossible];
	self.eigenfrequencies = [[lambda abs] sqrt]; self.eigenfrequencies.name = @"eigenfrequencies";
	[self normalizeEigenvectors:S withNorm: norm];
}

- (void) normalizeEigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm
{
    GLLinearTransform *T, *Tinv, *T_zIn;
    if (self.toDimensions) {
        T = [self.T expandedWithFromDimensions: self.fromDimensions toDimensions: self.toDimensions];
        T_zIn = [self.T_zIn expandedWithFromDimensions: self.fromDimensions toDimensions: self.toFinalDimensions];
        Tinv = [GLLinearTransform discreteTransformFromDimension: self.zDim toDimension: self.chebDim forEquation:self.equation];
        Tinv = [Tinv expandedWithFromDimensions: self.toDimensions toDimensions: self.fromDimensions];
    } else {
        T = self.T;
        T_zIn = self.T_zIn;
        Tinv = [GLLinearTransform discreteTransformFromDimension: self.zDim toDimension: self.chebDim forEquation:self.equation];
    }
    
    GLLinearTransform *S_spatial_cheb_grid = [T matrixMultiply: [S makeRealIfPossible]];
    S_spatial_cheb_grid = [S_spatial_cheb_grid normalizeWithFunction: norm];
    S = [Tinv matrixMultiply: S_spatial_cheb_grid];
    
	self.S = [T_zIn matrixMultiply: S];
	self.S.name = @"S_transform";
	
	GLLinearTransform *diffZ;
	if (S.toDimensions.count == diffZ.fromDimensions.count) {
		diffZ = self.diffZ;
	} else {
		diffZ = [self.diffZ expandedWithFromDimensions: S.toDimensions toDimensions:S.toDimensions];
	}
	
	self.Sprime = [T_zIn matrixMultiply: [diffZ multiply: S]];
	GLLinearTransform *scaling = [GLLinearTransform linearTransformFromFunction: self.eigendepths];
	GLMatrixMatrixDiagonalDenseMultiplicationOperation *op = [[GLMatrixMatrixDiagonalDenseMultiplicationOperation alloc] initWithFirstOperand: self.Sprime secondOperand: scaling];
	self.Sprime = op.result[0]; self.Sprime.name = @"Sprime_transform";
}

- (NSArray *) internalGeostrophicModesFromDensityProfile: (GLFunction *) rho forLatitude: (GLFloat) latitude
{
	[self createStratificationProfileFromDensity: rho atLatitude: latitude];
	
	GLLinearTransform *A = self.Tzz;
	GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [self.N2_cheb_grid times: @(-1/g)]] multiply: self.T];
	
	GLFloat *a = A.pointerValue;
	GLFloat *b = B.pointerValue;
	GLFloat *T = self.T.pointerValue;
	NSUInteger iTop=0;
	NSUInteger iBottom=A.matrixDescription.strides[0].nRows-1;
	for (NSUInteger j=0; j<A.matrixDescription.strides[0].nColumns; j++) {
		a[iTop*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = T[iTop*self.T.matrixDescription.strides[0].rowStride + j*self.T.matrixDescription.strides[0].columnStride];
		b[iTop*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = 0.0;
		
		a[iBottom*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = T[iBottom*self.T.matrixDescription.strides[0].rowStride + j*self.T.matrixDescription.strides[0].columnStride];
		b[iBottom*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = 0.0;
	}
	
	NSArray *system = [B generalizedEigensystemWith: A];
	
	[self normalizeDepthBasedEigenvalues: system[0] eigenvectors: system[1] withNorm: [self.N2_cheb_grid times: @(1/g)]];
	
	self.eigenfrequencies = [self.eigendepths times: @(0)];
	return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k forLatitude: (GLFloat) latitude
{
	[self createStratificationProfileFromDensity: rho atLatitude: latitude];
	
	GLScalar *k2 = [GLScalar scalarWithValue: k*k forEquation: self.equation];
	GLLinearTransform *A = [self.Tzz minus: [k2 multiply: self.T]];
    GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(-1/g)]] multiply: self.T];
	
    GLFloat *a = A.pointerValue;
    GLFloat *b = B.pointerValue;
    GLFloat *T = self.T.pointerValue;
    NSUInteger iTop=0;
    NSUInteger iBottom=A.matrixDescription.strides[0].nRows-1;
    for (NSUInteger j=0; j<A.matrixDescription.strides[0].nColumns; j++) {
        a[iTop*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = T[iTop*self.T.matrixDescription.strides[0].rowStride + j*self.T.matrixDescription.strides[0].columnStride];
        b[iTop*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = 0.0;
        
        a[iBottom*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = T[iBottom*self.T.matrixDescription.strides[0].rowStride + j*self.T.matrixDescription.strides[0].columnStride];
        b[iBottom*A.matrixDescription.strides[0].rowStride + j*A.matrixDescription.strides[0].columnStride] = 0.0;
    }
    
    NSArray *system = [B generalizedEigensystemWith: A];
	
	[self normalizeDepthBasedEigenvalues: system[0] eigenvectors: system[1] withNorm: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(1/g)]];
	
	self.eigenfrequencies = [[[[self.eigendepths abs] times: @(g*k*k)] plus: @(self.f0*self.f0)] sqrt];
	return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) createOptimizedOperationFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k forLatitude: (GLFloat) latitude maximumModes: (NSUInteger) maxModes zOutDim: (GLDimension *) zOutDim
{
    [self createStratificationProfileFromDensity: rho atLatitude: latitude];
    
    // These are the two matrices (without correct boundary conditions) that we need to find the eigensystem of.
    self.k2 = [GLScalar scalarWithValue: k*k forEquation: self.equation];
    GLLinearTransform *A = [self.Tzz minus: [self.k2 multiply: self.T]];
    GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(-1/g)]] multiply: self.T];
    
    // Operation that sets the boundary conditions for A.
    GLMatrixDescription *mdA = A.matrixDescription;
    GLMatrixDescription *mdT = self.T.matrixDescription;
    variableOperation boundaryAop = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        GLFloat *Ain = (GLFloat *) [operandArray[0] bytes];
        GLFloat *Tin = (GLFloat *) [operandArray[1] bytes];
        GLFloat *Aout = (GLFloat *) [resultArray[0] bytes];
        
        memcpy(Aout, Ain, mdA.nPoints*sizeof(GLFloat));
        
        NSUInteger iTop=0;
        NSUInteger iBottom=mdA.strides[0].nRows-1;
        for (NSUInteger j=0; j<mdA.strides[0].nColumns; j++) {
            Aout[iTop*mdA.strides[0].rowStride + j*mdA.strides[0].columnStride] = Tin[iTop*mdT.strides[0].rowStride + j*mdT.strides[0].columnStride];
            Aout[iBottom*mdA.strides[0].rowStride + j*mdA.strides[0].columnStride] = Tin[iBottom*mdT.strides[0].rowStride + j*mdT.strides[0].columnStride];
        }
    };
    GLVariableOperation *operationA = [[GLVariableOperation alloc] initWithResult: @[[GLVariable variableWithPrototype: A]] operand:@[A,self.T] buffers:@[] operation: boundaryAop];
    A = operationA.result[0];
    
    // Operation that sets the boundary conditions for B.
    GLMatrixDescription *mdB = B.matrixDescription;
    variableOperation boundaryBop = ^(NSArray *resultArray, NSArray *operandArray, NSArray *bufferArray) {
        GLFloat *Bin = (GLFloat *) [operandArray[0] bytes];
        GLFloat *Bout = (GLFloat *) [resultArray[0] bytes];
        
        memcpy(Bout, Bin, mdB.nPoints*sizeof(GLFloat));
        
        NSUInteger iTop=0;
        NSUInteger iBottom=mdB.strides[0].nRows-1;
        for (NSUInteger j=0; j<mdB.strides[0].nColumns; j++) {
            Bout[iTop*mdB.strides[0].rowStride + j*mdB.strides[0].columnStride] = 0.0;
            Bout[iBottom*mdB.strides[0].rowStride + j*mdB.strides[0].columnStride] = 0.0;
        }
    };
    GLVariableOperation *operationB = [[GLVariableOperation alloc] initWithResult: @[[GLVariable variableWithPrototype: B]] operand:@[B] buffers:@[] operation: boundaryBop];
    B = operationB.result[0];
    
    // Now solve the eigenvalue problem.
    NSArray *system = [B generalizedEigensystemWith: A];
    
    GLLinearTransform *S;
    GLFunction *h;
    if (maxModes) {
        h = [system[0] variableFromIndexRangeString: [NSString stringWithFormat: @"0:%lu", maxModes-1]];
        S = [system[1] reducedFromDimensions: [NSString stringWithFormat: @"0:%lu", maxModes-1] toDimension: @":"];
    } else {
        h = system[0];
        S = system[1];
    }
    
    // Hopefully these matrices came out real
    S = [S realPart];
    h = [h realPart];
    
    
    // Now convert to the spatial domain in order to normalize (note, we can get around this if needed!)
    GLLinearTransform *Tinv = [GLLinearTransform discreteTransformFromDimension: self.zDim toDimension: self.chebDim forEquation:self.equation];
    GLLinearTransform *S_spatial_cheb_grid = [self.T matrixMultiply: S];
    S_spatial_cheb_grid = [S_spatial_cheb_grid normalizeWithFunction: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(1/g)]];
    GLLinearTransform *S_cheb = [Tinv matrixMultiply: S_spatial_cheb_grid];
    
    GLLinearTransform *T_zOut;
    if (self.zInDim == zOutDim) {
        T_zOut = self.T_zIn;
    } else {
        T_zOut = [GLLinearTransform discreteTransformFromDimension: self.chebDim toDimension: zOutDim forEquation: self.equation];
    }
    
    S = [T_zOut matrixMultiply: S_cheb];
    
    GLLinearTransform *Sprime = [T_zOut matrixMultiply: [self.diffZ multiply: S_cheb]];
    GLLinearTransform *scaling = [GLLinearTransform linearTransformFromFunction: h];
    GLMatrixMatrixDiagonalDenseMultiplicationOperation *op = [[GLMatrixMatrixDiagonalDenseMultiplicationOperation alloc] initWithFirstOperand: Sprime secondOperand: scaling];
    Sprime = op.result[0];
    
    self.optimizer = [[GLOperationOptimizer alloc] initWithTopVariables: @[self.k2] bottomVariables: @[h, S, Sprime]];
    self.operation = self.optimizer.operationBlock;
    self.dataBuffers = [NSMutableArray array];
    for (GLBuffer *aBuffer in self.optimizer.internalDataBuffers) {
        [self.dataBuffers addObject: [[GLMemoryPool sharedMemoryPool] dataWithLength: aBuffer.numBytes]];
    }
	//GLOperationVisualizer *vis = [[GLOperationVisualizer alloc] initWithTopVariables: self.optimizer.topVariables. bottomVariables:<#(NSArray *)#>]
	
	//NSLog(optimizer)
    
    S.name = @"S_transform";
    Sprime.name = @"Sprime_transform";
    h.name = @"eigendepths";
    
    return @[h, S, Sprime];
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k forLatitude: (GLFloat) latitude maximumModes: (NSUInteger) maxModes zOutDim: (GLDimension *) zOutDim
{
	if (!self.operation) {
        NSArray *results = [self createOptimizedOperationFromDensityProfile: rho wavenumber: k forLatitude:latitude maximumModes:maxModes zOutDim:zOutDim];
        
        self.eigendepths = results[0];
        self.S = results[1];
        self.Sprime = results[2];
	}
	
	GLScalar *k2 = [GLScalar scalarWithValue: k*k forEquation: self.equation];
	self.operation( @[self.eigendepths.data, self.S.data,self.Sprime.data], @[k2.data], self.dataBuffers);
	
	self.rossbyRadius = [[[self.eigendepths times: @(g/(self.f0*self.f0))] abs] sqrt]; self.rossbyRadius.name = @"rossbyRadii";
	self.eigenfrequencies = [[[[self.eigendepths abs] times: @(g*k*k)] plus: @(self.f0*self.f0)] sqrt];
	return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude maximumModes: (NSUInteger) maxModes zOutDim: (GLDimension *) zOutDim
{
	if (dimensions.count != 3) {
		[NSException raise:@"InvalidDimensions" format: @"We are assuming exactly three dimensions"];
	}
	
    NSArray *results = [self createOptimizedOperationFromDensityProfile: rho wavenumber: 0.0 forLatitude:latitude maximumModes:maxModes zOutDim:zOutDim];
    GLFunction *eigendepths = results[0];
    GLLinearTransform *S = results[1];
    GLLinearTransform *Sprime = results[2];
    
	// Now we need to create two sets of dimensions:
	// 1. Our 'fromDimensions' will be modes for z, and exponentials for x,y
	// 2. Our 'toDimensions' will be physical z, and exponentials for x,y
	
	// create an array with the intended transformation (this is agnostic to dimension ordering).
	NSMutableArray *fromBasis = [NSMutableArray array];
	NSMutableArray *toBasis = [NSMutableArray array];
	NSMutableArray *fromSpatialDimensions = [NSMutableArray array];
	NSMutableArray *toSpatialDimensions = [NSMutableArray array];
	NSMutableArray *matrixFormat = [NSMutableArray array];
	NSUInteger xIndex = NSNotFound;
	NSUInteger yIndex = NSNotFound;
	NSUInteger zIndex = NSNotFound;
	for (GLDimension *dim in dimensions) {
		if (self.zInDim == dim) {
			[matrixFormat addObject: @(kGLDenseMatrixFormat)];
			
			[toBasis addObject: @([S.toDimensions[0] basisFunction])];
			[toSpatialDimensions addObject: S.toDimensions[0]];

			[fromBasis addObject: @([S.fromDimensions[0] basisFunction])];
			[fromSpatialDimensions addObject: S.fromDimensions[0]];
			zIndex = [dimensions indexOfObject: dim];
		} else {
			[matrixFormat addObject: @(kGLDiagonalMatrixFormat)];
			
			[toBasis addObject: @(kGLExponentialBasis)];
			[toSpatialDimensions addObject: dim];
			
			[fromBasis addObject: @(kGLExponentialBasis)];
			[fromSpatialDimensions addObject: dim];
			if (xIndex == NSNotFound) {
				xIndex = [dimensions indexOfObject: dim];
			} else {
				yIndex = [dimensions indexOfObject: dim];
			}
		}
	}
	
	self.fromDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: fromSpatialDimensions transformedToBasis: fromBasis];
	self.toDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: toSpatialDimensions transformedToBasis: toBasis];

	GLDimension *kDim, *lDim;
	for (GLDimension *dim in self.toDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			kDim = dim;
		} else if ( [dim.name isEqualToString: @"l"]) {
			lDim = dim;
		}
	}
	
	self.k = [[GLFunction functionOfRealTypeFromDimension: kDim withDimensions: self.toDimensions forEquation: self.equation] scalarMultiply: 2*M_PI];
	self.l = [[GLFunction functionOfRealTypeFromDimension: lDim withDimensions: self.toDimensions forEquation: self.equation] scalarMultiply: 2*M_PI];
	
	GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: kDim withDimensions: self.fromDimensions forEquation: self.equation] scalarMultiply: 2*M_PI];
	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: lDim withDimensions: self.fromDimensions forEquation: self.equation] scalarMultiply: 2*M_PI];
	GLFunction *K2 = [[k multiply: k] plus: [l multiply: l]];
	
	self.S = [GLLinearTransform transformOfType:kGLRealDataFormat withFromDimensions:self.fromDimensions toDimensions:self.toDimensions inFormat:matrixFormat forEquation:self.equation matrix: nil];
	self.Sprime = [GLLinearTransform transformOfType:kGLRealDataFormat withFromDimensions:self.fromDimensions toDimensions:self.toDimensions inFormat:matrixFormat forEquation:self.equation matrix: nil];
	self.eigendepths = [GLFunction functionOfRealTypeWithDimensions: self.fromDimensions forEquation:self.equation];
	
	GLMatrixDescription *mdK2 = K2.matrixDescription;
	GLMatrixDescription *mdS = self.S.matrixDescription;
	for (NSUInteger i=0; i<mdK2.strides[xIndex].nPoints; i++ ) {
		for (NSUInteger j=0; j<mdK2.strides[yIndex].nPoints; j++ ) {
			NSUInteger pos = i*mdK2.strides[xIndex].stride + j*mdK2.strides[yIndex].stride;
			self.k2.pointerValue[0] = K2.pointerValue[pos];
			
			self.operation( @[eigendepths.data, S.data, Sprime.data], @[self.k2.data], self.dataBuffers);
			
			// Copy the eigenvalues over
			const GLFloat *A = eigendepths.data.bytes;
			GLFloat *B = self.eigendepths.data.mutableBytes;
			NSUInteger offset =i*mdK2.strides[xIndex].stride + j*mdK2.strides[yIndex].stride;
			for (NSUInteger n=0; n<mdK2.strides[zIndex].nPoints; n++ ) {
				B[offset + n*mdK2.strides[zIndex].stride] = A[n];
			}
			
			const GLFloat *S1 = S.data.bytes;
			const GLFloat *Sprime1 = Sprime.data.bytes;
			GLFloat *S2 = self.S.data.mutableBytes;
			GLFloat *Sprime2 = self.Sprime.data.mutableBytes;
			
			offset =i*mdS.strides[xIndex].stride + j*mdS.strides[yIndex].stride;
			for (NSUInteger iPoint=0; iPoint<mdS.strides[zIndex].nPoints; iPoint++ ) {
				S2[offset + iPoint*mdS.strides[zIndex].stride] = S1[iPoint];
				Sprime2[offset + iPoint*mdS.strides[zIndex].stride] = Sprime1[iPoint];
			}
		}
	}
	
	
	self.eigenfrequencies = [[[[self.eigendepths abs] multiply: [K2 times: @(g)]] plus: @(self.f0*self.f0)] sqrt];
	
	return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude
{
    [self createStratificationProfileFromDensity: rho atLatitude: latitude];
    
    if (dimensions.count != 3) {
        [NSException raise:@"InvalidDimensions" format: @"We are assuming exactly three dimensions"];
    }
    
    // Now we need to create two sets of dimensions:
    // 1. Our 'fromDimensions' will be Chebyshev polynomails for z, and exponentials for x,y
    // 2. Our 'toDimensions' will be physical z, and exponentials for x,y
    
	// create an array with the intended transformation (this is agnostic to dimension ordering).
	NSMutableArray *fromBasis = [NSMutableArray array];
    NSMutableArray *toBasis = [NSMutableArray array];
    NSMutableArray *fromSpatialDimensions = [NSMutableArray array];
    NSMutableArray *toSpatialDimensions = [NSMutableArray array];
    NSMutableArray *toFinalSpatialDimensions = [NSMutableArray array];
    NSUInteger xIndex = NSNotFound;
    NSUInteger yIndex = NSNotFound;
    NSUInteger zIndex = NSNotFound;
	for (GLDimension *dim in dimensions) {
        if (self.zInDim == dim) {
            [toBasis addObject: @(self.zDim.basisFunction)];
            [toSpatialDimensions addObject: self.zDim];
            [toFinalSpatialDimensions addObject: self.zInDim];
            [fromBasis addObject: @(self.chebDim.basisFunction)];
            [fromSpatialDimensions addObject: self.chebDim];
            zIndex = [dimensions indexOfObject: dim];
        } else {
            [toBasis addObject: @(kGLExponentialBasis)];
            [toSpatialDimensions addObject: dim];
            [toFinalSpatialDimensions addObject: dim];
            [fromBasis addObject: @(kGLExponentialBasis)];
            [fromSpatialDimensions addObject: dim];
            if (xIndex == NSNotFound) {
                xIndex = [dimensions indexOfObject: dim];
            } else {
                yIndex = [dimensions indexOfObject: dim];
            }
        }
	}
	
	self.fromDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: fromSpatialDimensions transformedToBasis: fromBasis];
    self.toDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: toSpatialDimensions transformedToBasis: toBasis];
    self.toFinalDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: toFinalSpatialDimensions transformedToBasis: toBasis];
	GLDimension *kDim, *lDim;
	for (GLDimension *dim in self.toDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			kDim = dim;
		} else if ( [dim.name isEqualToString: @"l"]) {
			lDim = dim;
		}
	}
	
	GLEquation *equation = rho.equation;
	self.k = [[GLFunction functionOfRealTypeFromDimension: kDim withDimensions: self.toDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	self.l = [[GLFunction functionOfRealTypeFromDimension: lDim withDimensions: self.toDimensions forEquation: equation] scalarMultiply: 2*M_PI];
    GLLinearTransform *K2 = [[GLLinearTransform linearTransformFromFunction:[[self.k multiply: self.k] plus: [self.l multiply: self.l]]] expandedWithFromDimensions: self.toDimensions toDimensions: self.toDimensions];
	
	GLLinearTransform *T = [self.T expandedWithFromDimensions: self.fromDimensions toDimensions: self.toDimensions];
    GLLinearTransform *Tzz = [self.Tzz expandedWithFromDimensions: self.fromDimensions toDimensions: self.toDimensions];
    
    GLLinearTransform *A = [Tzz minus: [K2 multiply: T]];
    GLLinearTransform *scaling = [[GLLinearTransform linearTransformFromFunction: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(-1/g)]] expandedWithFromDimensions: self.toDimensions toDimensions: self.toDimensions];
    GLLinearTransform *B = [scaling multiply: T];
    
    GLFloat *a = A.pointerValue;
    GLFloat *b = B.pointerValue;
    GLMatrixDescription *md = A.matrixDescription;
    NSUInteger iTop=0;
    NSUInteger iBottom = md.strides[zIndex].nRows-1;
    for (NSUInteger iX=0; iX<md.strides[xIndex].nDiagonalPoints; iX++) {
        NSUInteger index = iX*md.strides[xIndex].stride;
        for (NSUInteger iY=0; iY<md.strides[yIndex].nDiagonalPoints; iY++) {
            NSUInteger index2 = index + iY*md.strides[yIndex].stride;
            for (NSUInteger j=0; j<md.strides[zIndex].nColumns; j++) {
                NSUInteger topIndex = index2 + iTop*md.strides[zIndex].rowStride + j*md.strides[zIndex].columnStride;
                a[topIndex] = self.T.pointerValue[iTop*self.T.matrixDescription.strides[0].rowStride + j*self.T.matrixDescription.strides[0].columnStride];
                b[iTop*B.matrixDescription.strides[zIndex].rowStride + j*B.matrixDescription.strides[zIndex].columnStride] = 0.0;
                
                NSUInteger bottomIndex = index2 + iBottom*md.strides[zIndex].rowStride + j*md.strides[zIndex].columnStride;
                a[bottomIndex] = self.T.pointerValue[iBottom*self.T.matrixDescription.strides[0].rowStride + j*self.T.matrixDescription.strides[0].columnStride];
                b[iBottom*B.matrixDescription.strides[zIndex].rowStride + j*B.matrixDescription.strides[zIndex].columnStride] = 0.0;
            }
        }
    }
    
    NSArray *system = [B generalizedEigensystemWith: A];
	
	[self normalizeDepthBasedEigenvalues: system[0] eigenvectors: system[1] withNorm: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(1/g)]];
	
	NSArray *spectralDimensions = self.eigendepths.dimensions;
	GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: kDim withDimensions: spectralDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: lDim withDimensions: spectralDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	GLFunction *K2_spectral = [[k multiply: k] plus: [l multiply: l]];
	self.eigenfrequencies = [[[[self.eigendepths abs] multiply: [K2_spectral times: @(g)]] plus: @(self.f0*self.f0)] sqrt];
	
	return @[self.eigendepths, self.S, self.Sprime];
}

@end
