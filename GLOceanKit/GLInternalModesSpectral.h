//
//  GLInternalModesSpectral.h
//  GLOceanKit
//
//  Created by Jeffrey J. Early on 3/5/15.
//  Copyright (c) 2015 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface GLInternalModesSpectral : NSObject

/** Compute the internal geostrophic (omega=0) modes.
 @param rho Density profile given as a function of z only.
 @param latitude Latitude at which the modes will be used. The only affects the Rossby radii.
 @returns A GLInternalModes object with N2, eigendepths, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalGeostrophicModesFromDensityProfile: (GLFunction *) rho forLatitude: (GLFloat) latitude;

/** Compute the internal wave modes for a given wavenumber.
 @param rho Density profile given as a function of z only.
 @param k Wavenumber given in radians/meter.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @returns A GLInternalModes object with N2, eigendepths, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k forLatitude: (GLFloat) latitude;

/** Compute the internal wave modes a set of wavenumbers.
 @param rho Density profile given as a function of z only.
 @param dimensions A set of vertical and horizontal dimensions. The same vertical dimensions (z) must be included, and at least one horizontal dimension.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @returns A GLInternalModes object with N2, eigendepths, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalWaveModesFromDensityProfile: (GLFunction *) rho withFullDimensions: (NSArray *) dimensions forLatitude: (GLFloat) latitude;

// Coriolis frequency (in radians!) given the latitude.
@property GLFloat f0;

/// Density profile
@property(strong) GLFunction *rho;

/// Stratification profile as computed for the mode calculation. N^2(z) = -g/mean(rho) * d/dz(rho)
@property(strong) GLFunction *N2;

/// Wavenumber function associated with x (in radians!). Horizontal dimensions are spectral, vertical is z. This may be nil.
@property(strong) GLFunction *k;

/// Wavenumber function associated with y (in radians!). Horizontal dimensions are spectral, vertical is z. This may be nil.
@property(strong) GLFunction *l;

/// The eigenvalue omega (will be zero for geostrophic modes)
@property(strong) GLFunction *eigenfrequencies;

/// The eigenvalue h, the equivalent depth.
@property(strong) GLFunction *eigendepths;

/// Rossby radii associated with each mode (given its equivalent depth and latitude)
@property(strong) GLFunction *rossbyRadius;

/// Transformation from the eigenbasis for the w-modes to z.
@property(strong) GLLinearTransform *S;

/// Transformation from the eigenbasis for the (u,v)-modes to z.
@property(strong) GLLinearTransform *Sprime;

@property(strong) GLLinearTransform *diffOp;

@property(readonly,strong) GLEquation *equation;

/** Compute the internal wave modes for a given wavenumber.
 @discussion This method differs from the above, primary methods, by solving the generalized eigenvalue from for omega, rather than h.
 @param rho Density profile given as a function of z only.
 @param k Wavenumber given in radians/meter.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @returns A GLInternalModes object with N2, eigendepths, eigenfrequencies, S, and Sprime populated.
 */
- (NSArray *) internalWaveModesUsingGEPFromDensityProfile: (GLFunction *) rho wavenumber: (GLFloat) k forLatitude: (GLFloat) latitude;

@end
