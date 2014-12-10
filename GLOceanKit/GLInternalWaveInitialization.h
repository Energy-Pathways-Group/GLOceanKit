//
//  GLInternalWaveInitialization.h
//  InternalWaves
//
//  Created by Jeffrey J. Early on 1/14/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import "GLInternalModes.h"

@interface GLInternalWaveInitialization : NSObject <NSCoding>

/** Compute the internal wave modes a set of wavenumbers.
 @param rho Density profile given as a function of z only.
 @param dimensions A set of vertical and horizontal dimensions. The same vertical dimensions (z) must be included, and at least one horizontal dimension.
 @param latitude Latitude at which the modes will be used (for the Coriolis frequency).
 @param equation The GLEquation object that sould be used for all calculations.
 @returns A GLInternalWaveInitialization object with N2, eigenfrequencies, S, and Sprime and all the variable phases populated.
 */
- (GLInternalWaveInitialization *) initWithDensityProfile: (GLFunction *) rho fullDimensions: (NSArray *) dimensions latitude: (GLFloat) latitude equation: (GLEquation *) equation;

/// Optional maximum number of modes to be included in the transformation matrices. Setting this to 0 will use all modes.
@property NSUInteger maximumModes;

// Optionally restrict the depth of the modes (after they've been compute for their full depth)
@property GLFloat minDepth;

// Optionally restrict the depth of the modes (after they've been compute for their full depth)
@property GLFloat maxDepth;

/** Initializes all variables with the Garrett-Munk spectrum.
 @param energyLevel A multiplicative factor, use 1.0 for default settings.
 */
- (void) createGarrettMunkSpectrumWithEnergy: (GLFloat) energyLevel;

// Returns omega, the frequency of the wave.
- (GLFloat) createUnitWaveWithSpeed: (GLFloat) U_max verticalMode: (NSUInteger) mode k: (NSUInteger) kUnit l: (NSUInteger) lUnit omegaSign: (GLFloat) sign;

/// The equation used for all computations.
@property(strong) GLEquation *equation;

/// The object used to generate the internal modes. This contains the untruncated (not limited to maximumModes) matrices.
@property(strong) GLInternalModes *internalModes;

/// The spatial dimensions, e.g., (x, y, z), although they will be in the order given during initialization.
@property(strong) NSArray *fullDimensions;

/// The associated spectral dimensions, e.g., (k, l, mode).
@property(strong) NSArray *spectralDimensions;

/// The latitude used for the mode and phase computation.
@property GLFloat latitude;

/// The Coriolis frequency at the given latitude.
@property GLFloat f0;

/// Density profile (function of z only).
@property(strong) GLFunction *rho;

/// Stratification profile used in the calculations.
@property(strong) GLFunction *N2;

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

@property(strong) GLFunction *zeta_plus;
@property(strong) GLFunction *zeta_minus;
@property(strong) GLFunction *rho_plus;
@property(strong) GLFunction *rho_minus;
@property(strong) GLFunction *u_plus;
@property(strong) GLFunction *u_minus;
@property(strong) GLFunction *v_plus;
@property(strong) GLFunction *v_minus;
@property(strong) GLFunction *w_plus;
@property(strong) GLFunction *w_minus;

@end
