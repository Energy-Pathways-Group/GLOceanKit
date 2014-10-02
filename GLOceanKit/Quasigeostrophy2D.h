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

// Current SSH. Either 0 if freshly initialized, or the the last recorded SSH value if initialized from file.
@property(strong) GLFunction *ssh;

/************************************************/
/*		Fixed Parameters						*/
/************************************************/

#pragma mark -
#pragma mark Fixed Parameters
#pragma mark

/// The equivalent depth
@property(readonly) GLFloat h;

/// The latitude used for the mode and phase computation.
@property(readonly) GLFloat latitude;

/// Coriolis frequency (in radians!) given the latitude.
@property(readonly) GLFloat f0;

/// Variation in the Coriolis parameter
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
@property BOOL shouldUseBeta;

// Whether or not to use spectral vanishing viscosity
@property BOOL shouldUseSVV;

// Whether or not to use the 2/3 anti-aliasing rules
@property BOOL shouldAntiAlias;

// What wavenumber do you want to force at, relative to...
@property GLFloat forcingFraction;

@property GLFloat thermalDampingFraction;
@property GLFloat frictionalDampingFraction;
@property GLFloat forcingWidth;
@property GLFloat f_zeta; // beta = 0.1
@property GLFloat forcingDecorrelationTime;

@end
