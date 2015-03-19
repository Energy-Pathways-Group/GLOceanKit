//
//  GLOceanKitTests.m
//  GLOceanKitTests
//
//  Created by Jeffrey J. Early on 10/2/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <XCTest/XCTest.h>
@import GLNumericalModelingKit;
@import GLOceanKit;
//#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
//#import <GLOceanKit/GLOceanKit.h>

@interface GLOceanKitTests : XCTestCase

@end

// First check if their difference is less than precision, then do the same, but scaled by the magnitude.
#define fequal(a,b) ((fabs((a) - (b)) < 10*FLT_EPSILON) || (fabs(((a) - (b))/a) < 10*FLT_EPSILON))

// Set your own precision
#define fequalprec(a,b,prec) ((fabs((a) - (b)) < prec) || (fabs(((a) - (b))/a) < prec))

@implementation GLOceanKitTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (void)testInternalModesFiniteDifference {
	GLFloat N2_0 = 1.69e-4;
	GLFloat rho0 = 1025;
	GLFloat g = 9.81;
    GLFloat latitude = 33.0;
    GLFloat H = 300;
	
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 32 domainMin: -H length: H];
	GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
	GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
	
	GLInternalModes *internalModes = [[GLInternalModes alloc] init];
	[internalModes internalGeostrophicModesFromDensityProfile: rho_bar forLatitude: latitude];

    GLLinearTransform *S = internalModes.S;
    GLLinearTransform *Sprime = internalModes.Sprime;
    
    GLDimension *nDim = S.fromDimensions[0];
    GLFunction *n = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim,zDim] forEquation: equation];
    GLFunction *z2 = [GLFunction functionOfRealTypeFromDimension: zDim withDimensions:@[nDim,zDim] forEquation: equation];
    
    GLFloat coeff = sqrt(2*g/H)/sqrt(N2_0-internalModes.f0*internalModes.f0);
    GLFunction *kz = [[n plus: @(1.)] times: @(M_PI/H)];
    GLFunction *hn = [[[[n plus: @(1.)] times: @(M_PI/H)] pow: -2.0] times: @(N2_0/g)];
    GLFunction *Gn = [[[kz times: z2] sin] times: @(coeff)];
    GLFunction *Fn = [[[[[kz times: z2] cos] times: @(coeff)] times: kz] times: hn];
    
    NSUInteger A_rs = Gn.matrixDescription.strides[1].stride;
    NSUInteger A_cs = Gn.matrixDescription.strides[0].stride;
    
    NSUInteger B_rs = S.matrixDescription.strides[0].rowStride;
    NSUInteger B_cs = S.matrixDescription.strides[0].columnStride;
    
    GLFloat precision = 1e-2;
    
    // First check the w-modes, Gn
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        GLFloat max = 0.0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if (fabs(Gn.pointerValue[i*A_rs + j*A_cs]) > max) {
                max = fabs(Gn.pointerValue[i*A_rs + j*A_cs]);
            }
        }
        
        NSUInteger failures = 0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if ( !fequalprec(fabs(Gn.pointerValue[i*A_rs + j*A_cs])/max, fabs(S.pointerValue[i*B_rs + j*B_cs])/max,precision) ) {
                //XCTFail(@"(%lu,%lu): Expected %f, found %f.",i,j, Gn.pointerValue[i*A_rs + j*A_cs], S.pointerValue[i*B_rs + j*B_cs]);
                failures++;
            }
        }
        if (failures) {
            XCTFail(@"w-mode %lu failed for %lu (of %lu) values with a precision of %f", j, failures, nDim.nPoints, precision);
            break;
        }
    }
    
    // Next check the (u,v)-modes, Fn
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        GLFloat max = 0.0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if (fabs(Fn.pointerValue[i*A_rs + j*A_cs]) > max) {
                max = fabs(Fn.pointerValue[i*A_rs + j*A_cs]);
            }
        }
        
        NSUInteger failures = 0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if ( !fequalprec(fabs(Fn.pointerValue[i*A_rs + j*A_cs])/max, fabs(Sprime.pointerValue[i*B_rs + j*B_cs])/max,precision) ) {
                //XCTFail(@"(%lu,%lu): Expected %f, found %f.",i,j, Fn.pointerValue[i*A_rs + j*A_cs], S.pointerValue[i*B_rs + j*B_cs]);
                failures++;
            }
        }
        if (failures) {
            XCTFail(@"(u,v)-mode %lu failed for %lu (of %lu) values with a precision of %f", j, failures, nDim.nPoints, precision);
            break;
        }
    }
    
    n = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim] forEquation: equation];
    GLFunction *h = [[[[n plus: @(1.)] times: @(M_PI/H)] pow: -2.0] times: @(N2_0/g)];
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        if (!fequalprec(fabs(internalModes.eigendepths.pointerValue[j])/h.pointerValue[j], 1.0,precision)) {
            XCTFail(@"Eigenvalue %lu failed with a precision of %f. Expected %g, found %g", j, precision, h.pointerValue[j], internalModes.eigendepths.pointerValue[j]);
            break;
        }
        
    }
    
}

- (void)testInternalModesSpectral {
    GLFloat N2_0 = 1.69e-4;
    GLFloat rho0 = 1025;
    GLFloat g = 9.81;
    GLFloat latitude = 33.0;
    GLFloat H = 300;
    
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 32 domainMin: -H length: H]; zDim.name = @"z";
    GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
    GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
    
    GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
    [internalModes internalGeostrophicModesFromDensityProfile: rho_bar forLatitude: latitude];
    
    GLLinearTransform *S = internalModes.S;
    GLLinearTransform *Sprime = internalModes.Sprime;
    
    GLDimension *nDim = S.fromDimensions[0];
    GLFunction *n = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim,zDim] forEquation: equation];
    GLFunction *z2 = [GLFunction functionOfRealTypeFromDimension: zDim withDimensions:@[nDim,zDim] forEquation: equation];
    
    GLFloat coeff = sqrt(2*g/H)/sqrt(N2_0-internalModes.f0*internalModes.f0);
    GLFunction *kz = [[n plus: @(1.)] times: @(M_PI/H)];
    GLFunction *hn = [[[[n plus: @(1.)] times: @(M_PI/H)] pow: -2.0] times: @(N2_0/g)];
    GLFunction *Gn = [[[kz times: z2] sin] times: @(coeff)];
    GLFunction *Fn = [[[[[kz times: z2] cos] times: @(coeff)] times: kz] times: hn];
    
    NSUInteger A_rs = Gn.matrixDescription.strides[1].stride;
    NSUInteger A_cs = Gn.matrixDescription.strides[0].stride;
    
    NSUInteger B_rs = S.matrixDescription.strides[0].rowStride;
    NSUInteger B_cs = S.matrixDescription.strides[0].columnStride;
    
    GLFloat precision = 1e-2;
    
    // First check the w-modes, Gn
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        GLFloat max = 0.0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if (fabs(Gn.pointerValue[i*A_rs + j*A_cs]) > max) {
                max = fabs(Gn.pointerValue[i*A_rs + j*A_cs]);
            }
        }
        
        NSUInteger failures = 0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if ( !fequalprec(fabs(Gn.pointerValue[i*A_rs + j*A_cs])/max, fabs(S.pointerValue[i*B_rs + j*B_cs])/max,precision) ) {
                //XCTFail(@"(%lu,%lu): Expected %f, found %f.",i,j, Gn.pointerValue[i*A_rs + j*A_cs], S.pointerValue[i*B_rs + j*B_cs]);
                failures++;
            }
        }
        if (failures) {
            XCTFail(@"w-mode %lu failed for %lu (of %lu) values with a precision of %f", j, failures, nDim.nPoints, precision);
            break;
        }
    }
    
    // Next check the (u,v)-modes, Fn
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        GLFloat max = 0.0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if (fabs(Fn.pointerValue[i*A_rs + j*A_cs]) > max) {
                max = fabs(Fn.pointerValue[i*A_rs + j*A_cs]);
            }
        }
        
        NSUInteger failures = 0;
        for (NSUInteger i=0; i<zDim.nPoints; i++) {
            if ( !fequalprec(fabs(Fn.pointerValue[i*A_rs + j*A_cs])/max, fabs(Sprime.pointerValue[i*B_rs + j*B_cs])/max,precision) ) {
                //XCTFail(@"(%lu,%lu): Expected %f, found %f.",i,j, Fn.pointerValue[i*A_rs + j*A_cs], S.pointerValue[i*B_rs + j*B_cs]);
                failures++;
            }
        }
        if (failures) {
            XCTFail(@"(u,v)-mode %lu failed for %lu (of %lu) values with a precision of %f", j, failures, nDim.nPoints, precision);
            break;
        }
    }
    
    n = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim] forEquation: equation];
    GLFunction *h = [[[[n plus: @(1.)] times: @(M_PI/H)] pow: -2.0] times: @(N2_0/g)];
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        if (!fequalprec(fabs(internalModes.eigendepths.pointerValue[j])/h.pointerValue[j], 1.0,precision)) {
            XCTFail(@"Eigenvalue %lu failed with a precision of %f. Expected %g, found %g", j, precision, h.pointerValue[j], internalModes.eigendepths.pointerValue[j]);
            break;
        }
        
    }
}

//- (void)testPerformanceExample {
//    // This is an example of a performance test case.
//    [self measureBlock:^{
//        // Put the code you want to measure the time of here.
//    }];
//}

@end
