//
//  Quasigeostrophy2D.h
//  GLOceanKit
//
//  Created by Jeffrey J. Early on 10/2/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Foundation/Foundation.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>

@interface Quasigeostrophy2D : NSObject

- (Quasigeostrophy2D *) initWithDimensions: (NSArray *) dims depth: (GLFloat) h latitude: (GLFloat) lat0 equation: (GLEquation *) equation;

- (Quasigeostrophy2D *) initWithFile: (NSURL *) fileURL resolutionDoubling: (BOOL) shouldDouble;


/// The spatial dimensions, e.g., (x, y), although they will be in the order given during initialization.
@property(strong) NSArray *dimensions;

@property(strong) GLEquation *equation;

@property(strong) GLRungeKuttaOperation *integrator;

// Current SSH. Either 0 everywhere, if freshly initialized, or the the last recorded SSH value, if initialized from file.
@property(strong) GLFunction *ssh;

// Initial phase of the forcing function
@property(strong) GLFunction *phi;

// Position of the floats
@property(strong) GLFunction *xPosition;
@property(strong) GLFunction *yPosition;

// An array of passive tracers
@property(strong) NSArray *tracers;

/************************************************/
/*		Fixed Parameters						*/
/************************************************/

#pragma mark -
#pragma mark Fixed Parameters
#pragma mark

/// The equivalent depth [meters].
@property(readonly) GLFloat h;

/// The latitude [degrees] used for the mode and phase computation.
@property(readonly) GLFloat latitude;

/// Coriolis frequency [radians] given the latitude.
@property(readonly) GLFloat f0;

/// Variation in the Coriolis parameter [radians/meter]
@property(readonly) GLFloat beta;

/// Length, time, and height scale
@property(readonly) GLFloat L_QG;
@property(readonly) GLFloat T_QG;
@property(readonly) GLFloat N_QG;

/************************************************/
/*		Adjustable Parameters					*/
/************************************************/

#pragma mark -
#pragma mark Adjustable Parameters
#pragma mark

// Whether or not to use beta-plane dynamics
// Defaults to NO.
@property BOOL shouldUseBeta;

// Whether or not to use spectral vanishing viscosity
// Defaults to YES.
@property BOOL shouldUseSVV;

// Whether or not to use the 2/3 anti-aliasing rules
// Defaults to NO.
@property(nonatomic) BOOL shouldAntiAlias;

// Whether or not to use narrow band forcing
// Defaults to YES.
@property BOOL shouldForce;

// What wavenumber do you want to force at, relative to the largest wavenumber.
// k_f = k_max/forcingFraction.
// Note that this value automatically doubles if you initialize with resolution doubling.
@property GLFloat forcingFraction;

// Width of the forcing, relative to the smallest wavenumber (wavenumber units).
@property GLFloat forcingWidth;

// Forcing strength, in nondimensional units.
@property GLFloat f_zeta; // beta = 0.1

// Forcing decorrelation time in days. HUGE_VAL means the phases don't change.
@property GLFloat forcingDecorrelationTime;

// What wavenumber do you want to damp at, relative to the smallest wavenumber.
// k_damp = k_min*dampingFraction
@property GLFloat thermalDampingFraction;

// What wavenumber do you want to damp at, relative to the smallest wavenumber.
// k_damp = k_min*dampingFraction
@property GLFloat frictionalDampingFraction;

@end
