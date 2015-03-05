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
@property(strong) GLDimension *zDim;
@property(strong) GLLinearTransform *diffZ;
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
		self.zDim = zInitial;
		self.rho = rho;
	} else {
		self.zDim = [[GLDimension alloc] initDimensionWithGrid: kGLChebyshevEndpointGrid nPoints: zInitial.nPoints domainMin:zInitial.domainMin length:zInitial.domainLength];
		GLFunction *z = [GLFunction functionOfRealTypeFromDimension: self.zDim withDimensions: @[self.zDim] forEquation: self.equation];
		self.rho = [rho interpolateAtPoints: @[z]];
	}
	
	GLFunction *rho_cheb = [self.rho transformToBasis: @[@(kGLChebyshevBasis)]];
	GLFunction *buoyancy = [[rho dividedBy: rho0] times: @(-g)];
	
	self.diffZ = [GLLinearTransform differentialOperatorWithDerivatives: 1 fromDimension: buoyancy.dimensions[0] forEquation: self.equation];
	
	// First construct N^2
	self.diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLNeumannBoundaryCondition rightBC:kGLNeumannBoundaryCondition bandwidth:1 fromDimension:self.zDim forEquation: self.equation];
	//self.diffZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 1 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension:self.zDim forEquation: self.equation];
	self.N2 = [self.diffZ transform: [[rho dividedBy: rho0] times: @(-g)]];
	self.N2.name = @"N2";
}

- (void) normalizeDepthBasedEigenvalues: (GLFunction *) lambda eigenvectors: (GLLinearTransform *) S withNorm: (GLFunction *) norm
{
	lambda = [lambda makeRealIfPossible];
	self.eigendepths = [lambda scalarDivide: 1.0]; self.eigendepths.name = @"eigendepths";
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
	S = [S makeRealIfPossible];
	self.S = [S normalizeWithFunction: norm]; self.S.name = @"S_transform";
	
	GLLinearTransform *diffZ;
	if (S.toDimensions.count == diffZ.fromDimensions.count) {
		diffZ = self.diffZ;
	} else {
		diffZ = [self.diffZ expandedWithFromDimensions: S.toDimensions toDimensions:S.toDimensions];
	}
	
	self.Sprime = [diffZ multiply: self.S];
	GLLinearTransform *scaling = [GLLinearTransform linearTransformFromFunction: self.eigendepths];
	GLMatrixMatrixDiagonalDenseMultiplicationOperation *op = [[GLMatrixMatrixDiagonalDenseMultiplicationOperation alloc] initWithFirstOperand: self.Sprime secondOperand: scaling];
	self.Sprime = op.result[0]; self.Sprime.name = @"Sprime_transform";
}

- (NSArray *) internalGeostrophicModesFromDensityProfile: (GLFunction *) rho forLatitude: (GLFloat) latitude
{
	[self createStratificationProfileFromDensity: rho atLatitude: latitude];
	
	GLFunction *invN2 = [self.N2 scalarDivide: -g];
	
	GLLinearTransform *diffZZ = [GLLinearTransform finiteDifferenceOperatorWithDerivatives: 2 leftBC: kGLDirichletBoundaryCondition rightBC:kGLDirichletBoundaryCondition bandwidth:1 fromDimension: self.zDim forEquation: self.equation];
	GLLinearTransform *invN2_trans = [GLLinearTransform linearTransformFromFunction: invN2];
	GLLinearTransform *diffOp = [invN2_trans multiply: diffZZ];
	
	NSArray *system = [diffOp eigensystemWithOrder: NSOrderedAscending];
	
	[self normalizeDepthBasedEigenvalues: system[0] eigenvectors: system[1] withNorm: [self.N2 times: @(1/g)]];
	
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
