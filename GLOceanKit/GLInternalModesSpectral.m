//
//  GLInternalModesSpectral.m
//  GLOceanKit
//
//  Created by Jeffrey J. Early on 3/5/15.
//  Copyright (c) 2015 Jeffrey J. Early. All rights reserved.
//

#import "GLInternalModesSpectral.h"
#import <GLNumericalModelingKit/GLLinearTransformationOperations.h>

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
    GLLinearTransform *S_spatial_cheb_grid = [self.T matrixMultiply: [S makeRealIfPossible]];
    S_spatial_cheb_grid = [S_spatial_cheb_grid normalizeWithFunction: norm];
    GLLinearTransform *Tinv = [GLLinearTransform discreteTransformFromDimension: self.zDim toDimension: self.chebDim forEquation:self.equation];
    S = [Tinv matrixMultiply: S_spatial_cheb_grid];
    
	self.S = [self.T_zIn matrixMultiply: S];
	self.S.name = @"S_transform";
	
	GLLinearTransform *diffZ;
	if (S.toDimensions.count == diffZ.fromDimensions.count) {
		diffZ = self.diffZ;
	} else {
		diffZ = [self.diffZ expandedWithFromDimensions: S.toDimensions toDimensions:S.toDimensions];
	}
	
	self.Sprime = [self.T_zIn matrixMultiply: [diffZ multiply: S]];
	GLLinearTransform *scaling = [GLLinearTransform linearTransformFromFunction: self.eigendepths];
	GLMatrixMatrixDiagonalDenseMultiplicationOperation *op = [[GLMatrixMatrixDiagonalDenseMultiplicationOperation alloc] initWithFirstOperand: self.Sprime secondOperand: scaling];
	self.Sprime = op.result[0]; self.Sprime.name = @"Sprime_transform";
}

- (NSArray *) internalGeostrophicModesFromDensityProfile: (GLFunction *) rho forLatitude: (GLFloat) latitude
{
	[self createStratificationProfileFromDensity: rho atLatitude: latitude];
	
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
	
	GLLinearTransform *A = self.Tzz;
	
//  GLLinearTransform *B = [[GLLinearTransform linearTransformFromFunction: [[self.N2_cheb_grid minus: @(self.f0*self.f0)] times: @(-g)]] multiply: self.T];
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
	
	// -g/(N2-f*f)
	GLFunction *invN2 = [[self.N2 minus: @(self.f0*self.f0)] scalarDivide: -g];
	
	GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:self.zDim forEquation:self.equation];
	GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
	GLLinearTransform *k2 = [GLLinearTransform transformOfType: kGLRealDataFormat withFromDimensions: @[self.zDim] toDimensions: @[self.zDim] inFormat: @[@(kGLDiagonalMatrixFormat)] forEquation:self.equation matrix:^( NSUInteger *row, NSUInteger *col ) {
		return (GLFloatComplex) (row[0]==col[0] ? k*k : 0);
	}];
	GLLinearTransform *diffOp = [invN2_trans multiply: [diffZZ minus: k2]];
	self.diffOp = diffOp;
	
	NSArray *system = [diffOp eigensystemWithOrder: NSOrderedAscending];
	
	[self normalizeDepthBasedEigenvalues: system[0] eigenvectors: system[1] withNorm: [[self.N2 minus: @(self.f0*self.f0)] times: @(1/g)]];
	
	self.eigenfrequencies = [[[[self.eigendepths abs] times: @(g*k*k)] plus: @(self.f0*self.f0)] sqrt];
	return @[self.eigendepths, self.S, self.Sprime];
}

- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude
{
	// create an array with the intended transformation (this is agnostic to dimension ordering).
	NSMutableArray *basis = [NSMutableArray array];
	GLDimension *zDim;
	for (GLDimension *dim in dimensions) {
		if ( [dim.name isEqualToString: @"x"] || [dim.name isEqualToString: @"y"]) {
			[basis addObject: @(kGLExponentialBasis)];
		} else {
			zDim = dim;
			[basis addObject: @(dim.basisFunction)];
		}
	}
	
	NSArray *transformedDimensions = [GLDimension dimensionsForRealFunctionWithDimensions: dimensions transformedToBasis: basis];
	GLDimension *kDim, *lDim;
	for (GLDimension *dim in transformedDimensions) {
		if ( [dim.name isEqualToString: @"k"]) {
			kDim = dim;
		} else if ( [dim.name isEqualToString: @"l"]) {
			lDim = dim;
		}
	}
	
	GLEquation *equation = rho.equation;
	self.k = [[GLFunction functionOfRealTypeFromDimension: kDim withDimensions: transformedDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	self.l = [[GLFunction functionOfRealTypeFromDimension: lDim withDimensions: transformedDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	GLFunction *K2 = [[self.k multiply: self.k] plus: [self.l multiply: self.l]];
	
	[self createStratificationProfileFromDensity: rho atLatitude: latitude];
	
	GLFunction *invN2 = [[self.N2 minus: @(self.f0*self.f0)] scalarDivide: -g];
	
	// Now construct A = k*k*eye(N) - Diff2;
	GLLinearTransform *diffZZ1D = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:zDim forEquation:equation];
	GLLinearTransform *diffZZ = [diffZZ1D expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];
	
	GLLinearTransform *invN2_trans = [[GLLinearTransform linearTransformFromFunction: invN2] expandedWithFromDimensions: transformedDimensions toDimensions: transformedDimensions];
	GLLinearTransform *diffOp = [invN2_trans multiply: [diffZZ minus: [GLLinearTransform linearTransformFromFunction:K2]]];
	
	NSArray *system = [diffOp eigensystemWithOrder: NSOrderedAscending];
	
	[self normalizeDepthBasedEigenvalues: system[0] eigenvectors: system[1] withNorm: [[self.N2 minus: @(self.f0*self.f0)] times: @(1/g)]];
	
	NSArray *spectralDimensions = self.eigendepths.dimensions;
	GLFunction *k = [[GLFunction functionOfRealTypeFromDimension: kDim withDimensions: spectralDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	GLFunction *l = [[GLFunction functionOfRealTypeFromDimension: lDim withDimensions: spectralDimensions forEquation: equation] scalarMultiply: 2*M_PI];
	GLFunction *K2_spectral = [[k multiply: k] plus: [l multiply: l]];
	self.eigenfrequencies = [[[[self.eigendepths abs] multiply: [K2_spectral times: @(g)]] plus: @(self.f0*self.f0)] sqrt];
	
	return @[self.eigendepths, self.S, self.Sprime];
}

@end
