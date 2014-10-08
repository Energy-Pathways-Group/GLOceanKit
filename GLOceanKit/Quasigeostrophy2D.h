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
@property(strong,readonly) NSArray *dimensions;
@property(strong,readonly) NSArray *wavenumberDimensions;
@property(strong,readonly) GLMutableDimension *tDim;

@property(strong,readonly) GLEquation *equation;

@property(strong,readonly) GLRungeKuttaOperation *integrator;

// These are automatically set for you, but you can override them.

// Current SSH. Either 0 everywhere, if freshly initialized, or the the last recorded SSH value, if initialized from file.
// Has dimensions that are given in the dimensions array.
@property(strong) GLFunction *ssh;

// Initial phase of the forcing function
// Has dimensions that are given in the wavenumberDimensions array.
@property(strong) GLFunction *phi;

// The variable 'forcing' contains the magnitude of the forcing term for each wavenumber, but no phase information.
// Has dimensions that are given in the wavenumberDimensions array.
@property(strong) GLFunction *forcing;

// the phase speed of each component
// Has dimensions that are given in the wavenumberDimensions array.
@property(strong) GLFunction *phaseSpeed;

// Position of the floats
// Can have any dimensions, but the values must refer the dimensions in the dimensions array.
@property(strong) GLFunction *xPosition;
@property(strong) GLFunction *yPosition;

// An array of passive tracers
// Has dimensions that are given in the dimensions array.
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

/// Whether or not you want to advect the floats specified by initial positions (xPosition, yPosition)
@property BOOL shouldAdvectFloats;

/// Whether or not you want to advect the tracers in the tracers array.
@property BOOL shouldAdvectTracer;

/************************************************/
/*		Output Control							*/
/************************************************/

#pragma mark -
#pragma mark Output Control
#pragma mark

// Adds the QG turbulence metadata to file. Do NOT call this function if you've set an outputFile below; that will be done automatically.
- (void) addMetadataToNetCDFFile: (GLNetCDFFile *) file;

/// Optional output file.
// By default this will output the ssh and any tracers or floats.
@property(copy) NSURL *outputFile;

/// Set to YES to write the ssh in frequency domain out. This is redundant and can be derived from the ssh.
@property BOOL shouldWriteSSHFD;

/// Set to YES to write the relative vorticity out. This is redundant and can be derived from the ssh.
@property BOOL shouldWriteRV;

/// Set to YES to write the force out. This is redundant and can be derived from the force magnitude and phase.
@property BOOL shouldWriteForce;

/// *Dimensionalized* ssh function (with time dimension) associated with the netcdf file.
@property(strong,readonly) GLMutableVariable *sshHistory;
@property(strong,readonly) GLFunction *dimensionalForceMag;
@property(strong,readonly) GLMutableVariable *phaseHistory;
@property(strong,readonly) GLMutableVariable *sshFDHistory;
@property(strong,readonly) GLMutableVariable *rvHistory;
@property(strong,readonly) GLMutableVariable *forceHistory;
@property(strong,readonly) GLMutableVariable *xPositionHistory;
@property(strong,readonly) GLMutableVariable *yPositionHistory;
@property(strong,readonly) NSArray *tracerHistories;

/************************************************/
/*		Derived Parameters						*/
/************************************************/

#pragma mark -
#pragma mark Derived Parameters
#pragma mark

// forcing wavenumber
@property(readonly) GLFloat k_f;

// small scale damping wavenumber
@property(readonly) GLFloat k_nu;

// large scale thermal damping wavenumber
@property(readonly) GLFloat k_alpha;

// large scale frictional damping wavenumber
@property(readonly) GLFloat k_r;

// forcing width
@property(readonly) GLFloat k_width;

// max resolved wavenumber
@property(readonly) GLFloat k_max;

// thermal damping parameter
@property(readonly) GLFloat alpha;

// frictional damping parameter
@property(readonly) GLFloat r;

// small scale damping parameter
@property(readonly) GLFloat nu;

@end
