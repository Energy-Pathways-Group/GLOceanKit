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

- (void) compareSolution: (GLInternalModes *) internalModes toAnalyticalSolutionWithK: (GLFloat) k N2: (GLFloat) N2_0 latitude: (GLFloat) latitude depth: (GLFloat) H
{
    GLFloat g = 9.81;
    
    GLLinearTransform *S = internalModes.S;
    GLLinearTransform *Sprime = internalModes.Sprime;
    
    GLDimension *zDim = S.toDimensions[0];
    GLDimension *nDim = S.fromDimensions[0];
    GLFunction *n = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim,zDim] forEquation: internalModes.equation];
    GLFunction *z2 = [GLFunction functionOfRealTypeFromDimension: zDim withDimensions:@[nDim,zDim] forEquation: internalModes.equation];
    
    GLFloat coeff = sqrt(2*g/H)/sqrt(N2_0-internalModes.f0*internalModes.f0);
    GLFunction *kz = [[n plus: @(1.)] times: @(M_PI/H)];
    GLFunction *hn = [[[[kz times: kz] plus: @(k*k)] pow: -1.0] times: @(N2_0/g)];
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
    
    n = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim] forEquation: internalModes.equation];
    kz = [[n plus: @(1.)] times: @(M_PI/H)];
    GLFunction *h = [[[[kz times: kz] plus: @(k*k)] pow: -1.0] times: @(N2_0/g)];
    for (NSUInteger j=0; j<nDim.nPoints; j++) {
        if (!fequalprec(fabs(internalModes.eigendepths.pointerValue[j])/h.pointerValue[j], 1.0,precision)) {
            XCTFail(@"Eigenvalue %lu failed with a precision of %f. Expected %g, found %g", j, precision, h.pointerValue[j], internalModes.eigendepths.pointerValue[j]);
            break;
        }
        
    }
}

- (void) compareSolution: (GLInternalModes *) internalModes toAnalyticalSolutionWithN2: (GLFloat) N2_0 latitude: (GLFloat) latitude depth: (GLFloat) H
{
    GLFloat g = 9.81;
    
    GLLinearTransform *S = internalModes.S;
    GLLinearTransform *Sprime = internalModes.Sprime;
    
    NSUInteger kIndex = NSNotFound, lIndex = NSNotFound, zIndex = NSNotFound;
    GLDimension *kDim, *lDim, *zDim, *nDim;
    for (NSUInteger i=0; i<S.matrixDescription.nDimensions; i++) {
        if (S.matrixDescription.strides[i].matrixFormat == kGLDenseMatrixFormat) {
            zIndex = i;
            zDim = S.toDimensions[zIndex];
            nDim = S.fromDimensions[zIndex];
        } else {
            if (!kDim) {
                kDim = S.toDimensions[i];
                kIndex = i;
            } else {
                lDim = S.toDimensions[i];
                lIndex = i;
            }
        }
    }
    GLFunction *kFunc = [GLFunction functionOfRealTypeFromDimension: kDim withDimensions:@[kDim] forEquation: internalModes.equation];
    GLFunction *lFunc = [GLFunction functionOfRealTypeFromDimension: lDim withDimensions:@[lDim] forEquation: internalModes.equation];
    GLFunction *zFunc = [GLFunction functionOfRealTypeFromDimension: zDim withDimensions:@[zDim] forEquation: internalModes.equation];
    GLFunction *nFunc = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:@[nDim] forEquation: internalModes.equation];
    
    GLFloat *k = (GLFloat *) kFunc.data.bytes;
    GLFloat *l = (GLFloat *) lFunc.data.bytes;
    GLFloat *z = (GLFloat *) zFunc.data.bytes;
    GLFloat *n = (GLFloat *) nFunc.data.bytes;
    
    GLFloat coeff = sqrt(2*g/H)/sqrt(N2_0-internalModes.f0*internalModes.f0);
    GLLinearTransform *Gn = [GLLinearTransform transformOfType: S.dataFormat withFromDimensions: S.fromDimensions toDimensions: S.toDimensions inFormat: S.matrixFormats forEquation: S.equation matrix:^( NSUInteger *row, NSUInteger *col ) {
        if ((row[kIndex] != col[kIndex]) || (row[lIndex] != col[lIndex]) ) return (GLFloatComplex) 0.0;
        
        GLFloat kz = (n[col[zIndex]] + 1.0)*M_PI/H;
        GLFloat zVal = z[row[zIndex]];
        
        GLFloatComplex value = coeff*sin(kz*zVal);
        
        return value;
    }];
    
    GLLinearTransform *Fn = [GLLinearTransform transformOfType: S.dataFormat withFromDimensions: S.fromDimensions toDimensions: S.toDimensions inFormat: S.matrixFormats forEquation: S.equation matrix:^( NSUInteger *row, NSUInteger *col ) {
        if ((row[kIndex] != col[kIndex]) || (row[lIndex] != col[lIndex]) ) return (GLFloatComplex) 0.0;
        
        GLFloat kz = (n[col[zIndex]] + 1.0)*M_PI/H;
        GLFloat k2 = (k[row[kIndex]]*k[row[kIndex]] + l[row[lIndex]]*l[row[lIndex]])*4.0*M_PI*M_PI;
        GLFloat h = (N2_0/g)/(kz*kz+k2);
        GLFloat zVal = z[row[zIndex]];
        
        GLFloatComplex value = coeff*kz*h*cos(kz*zVal);
        
        return value;
    }];
    
    GLFloat precision = 1e-2;
    for (NSUInteger iK=0; iK<S.matrixDescription.strides[kIndex].nDiagonalPoints; iK++) {
        NSUInteger index1 = iK*S.matrixDescription.strides[kIndex].stride;
        for (NSUInteger iL=0; iL<S.matrixDescription.strides[lIndex].nDiagonalPoints; iL++) {
            NSUInteger index2 = index1 + iL*S.matrixDescription.strides[lIndex].stride;
            for (NSUInteger iMode=0; iMode<S.matrixDescription.strides[zIndex].nColumns; iMode++) {
                NSUInteger index3 = index2 + iMode*S.matrixDescription.strides[zIndex].columnStride;
                
                GLFloat sModeMax = 0.0;
                for (NSUInteger iZ=0; iZ<S.matrixDescription.strides[zIndex].nRows; iZ++) {
                    NSUInteger index = index3 + iZ*S.matrixDescription.strides[zIndex].rowStride;
                    if (fabs(Gn.pointerValue[index]) > sModeMax) {
                        sModeMax = fabs(Gn.pointerValue[index]);
                    }
                }
                
                NSUInteger sModeFailures = 0;
                for (NSUInteger iZ=0; iZ<S.matrixDescription.strides[zIndex].nRows; iZ++) {
                    NSUInteger index = index3 + iZ*S.matrixDescription.strides[zIndex].rowStride;
                    if ( !fequalprec(fabs(Gn.pointerValue[index])/sModeMax, fabs(S.pointerValue[index])/sModeMax,precision) ) {
                        sModeFailures++;
                    }
                }
                
                if (sModeFailures) {
                    XCTFail(@"w-mode %lu @ (k,l)=(%lu,%lu) failed for %lu (of %lu) values with a precision of %f", iMode, iK, iL, sModeFailures, S.matrixDescription.strides[zIndex].nRows, precision);
                    break;
                }
            }
        }
    }
    
    for (NSUInteger iK=0; iK<S.matrixDescription.strides[kIndex].nDiagonalPoints; iK++) {
        NSUInteger index1 = iK*S.matrixDescription.strides[kIndex].stride;
        for (NSUInteger iL=0; iL<S.matrixDescription.strides[lIndex].nDiagonalPoints; iL++) {
            NSUInteger index2 = index1 + iL*S.matrixDescription.strides[lIndex].stride;
            for (NSUInteger iMode=0; iMode<S.matrixDescription.strides[zIndex].nColumns; iMode++) {
                NSUInteger index3 = index2 + iMode*S.matrixDescription.strides[zIndex].columnStride;
                
                GLFloat sModeMax = 0.0;
                for (NSUInteger iZ=0; iZ<S.matrixDescription.strides[zIndex].nRows; iZ++) {
                    NSUInteger index = index3 + iZ*S.matrixDescription.strides[zIndex].rowStride;
                    if (fabs(Fn.pointerValue[index]) > sModeMax) {
                        sModeMax = fabs(Fn.pointerValue[index]);
                    }
                }
                
                NSUInteger sModeFailures = 0;
                for (NSUInteger iZ=0; iZ<S.matrixDescription.strides[zIndex].nRows; iZ++) {
                    NSUInteger index = index3 + iZ*S.matrixDescription.strides[zIndex].rowStride;
                    if ( !fequalprec(fabs(Fn.pointerValue[index])/sModeMax, fabs(Sprime.pointerValue[index])/sModeMax,precision) ) {
                        sModeFailures++;
                    }
                }
                
                if (sModeFailures) {
                    XCTFail(@"(u,v)-mode %lu @ (k,l)=(%lu,%lu) failed for %lu (of %lu) values with a precision of %f", iMode, iK, iL, sModeFailures, S.matrixDescription.strides[zIndex].nRows, precision);
                    break;
                }
            }
        }
    }
    
    kFunc = [GLFunction functionOfRealTypeFromDimension: kDim withDimensions:internalModes.eigendepths.dimensions forEquation: internalModes.equation];
    lFunc = [GLFunction functionOfRealTypeFromDimension: lDim withDimensions:internalModes.eigendepths.dimensions forEquation: internalModes.equation];
    nFunc = [GLFunction functionOfRealTypeFromDimension: nDim withDimensions:internalModes.eigendepths.dimensions forEquation: internalModes.equation];
    
    GLFunction *kz = [[nFunc plus: @(1.)] times: @(M_PI/H)];
    GLFunction *K2 = [[[kFunc multiply: kFunc] plus: [lFunc multiply: lFunc]] times: @(M_PI*M_PI*4.0)];
    GLFunction *h = [[[[kz multiply: kz] plus: K2] pow: -1.0] times: @(N2_0/g)];
    
    for (NSUInteger iK=0; iK<h.matrixDescription.strides[kIndex].nPoints; iK++) {
        NSUInteger index1 = iK*h.matrixDescription.strides[kIndex].stride;
        for (NSUInteger iL=0; iL<h.matrixDescription.strides[lIndex].nPoints; iL++) {
            NSUInteger index2 = index1 + iL*h.matrixDescription.strides[lIndex].stride;
            for (NSUInteger iMode=0; iMode<h.matrixDescription.strides[zIndex].nPoints; iMode++) {
                NSUInteger index3 = index2 + iMode*h.matrixDescription.strides[zIndex].columnStride;
                if (!fequalprec(fabs(internalModes.eigendepths.pointerValue[index3])/h.pointerValue[index3], 1.0,precision)) {
                    XCTFail(@"Eigendepth %lu @ (k,l)=(%lu,%lu) failed with a precision of %f. Expected %g, found %g", iMode, iK, iL, precision, h.pointerValue[index3], internalModes.eigendepths.pointerValue[index3]);
                    break;
                }
            }
        }
    }
}

- (void)testInternalGeostrophicModesFiniteDifference {
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

    [self compareSolution: internalModes toAnalyticalSolutionWithK: 0.0 N2:N2_0 latitude:latitude depth:H];
    
}

- (void)testInternalGeostrophicModesSpectral {
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
    
    [self compareSolution: internalModes toAnalyticalSolutionWithK: 0.0 N2:N2_0 latitude:latitude depth:H];
}

- (void)testInternalWaveModesFiniteDifference {
    GLFloat N2_0 = 1.69e-4;
    GLFloat rho0 = 1025;
    GLFloat g = 9.81;
    GLFloat latitude = 33.0;
    GLFloat H = 300;
    GLFloat k = 0.1;
    
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 32 domainMin: -H length: H];
    GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
    GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
    
    GLInternalModes *internalModes = [[GLInternalModes alloc] init];
    [internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: k forLatitude: latitude];
    
    [self compareSolution: internalModes toAnalyticalSolutionWithK: k N2:N2_0 latitude:latitude depth:H];
    
}

- (void)testInternalWaveModesSpectral {
    GLFloat N2_0 = 1.69e-4;
    GLFloat rho0 = 1025;
    GLFloat g = 9.81;
    GLFloat latitude = 33.0;
    GLFloat H = 300;
    GLFloat k = 0.1;
    
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 512 domainMin: -H length: H]; zDim.name = @"z";
    GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
    GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
    
    GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
    [internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: k forLatitude: latitude];
    
    [self compareSolution: internalModes toAnalyticalSolutionWithK: k N2:N2_0 latitude:latitude depth:H];
}

- (void)testInternalWaveModesSpectralOptimized {
	GLFloat N2_0 = 1.69e-4;
	GLFloat rho0 = 1025;
	GLFloat g = 9.81;
	GLFloat latitude = 33.0;
	GLFloat H = 300;
	GLFloat k = 0.1;
	
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 512 domainMin: -H length: H]; zDim.name = @"z";
	GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
	GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
	
	GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
	[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: k forLatitude: latitude maximumModes: 31 zOutDim: zDim];
	
	[self compareSolution: internalModes toAnalyticalSolutionWithK: k N2:N2_0 latitude:latitude depth:H];
}

- (void)testInternalWaveFullModesFiniteDifference {
    GLFloat N2_0 = 1.69e-4;
    GLFloat rho0 = 1025;
    GLFloat g = 9.81;
    GLFloat latitude = 33.0;
    GLFloat H = 300;
    
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 32 domainMin: -H length: H]; zDim.name = @"z";
    GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
    GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
    
    GLFloat width = 1000;
    GLFloat height = 500;
    NSUInteger Nx = 4;
    NSUInteger Ny = 4;
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
    xDim.name = @"x";
    GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
    yDim.name = @"y";
    
    GLInternalModes *internalModes = [[GLInternalModes alloc] init];
    [internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude];
    
//    [internalModes.S dumpToConsole];
    
    [self compareSolution: internalModes toAnalyticalSolutionWithN2: N2_0 latitude: latitude depth: H];
    
}

- (void)testInternalWaveFullModesSpectral {
    GLFloat N2_0 = 1.69e-4;
    GLFloat rho0 = 1025;
    GLFloat g = 9.81;
    GLFloat latitude = 33.0;
    GLFloat H = 300;
    
    GLEquation *equation = [[GLEquation alloc] init];
    GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 32 domainMin: -H length: H]; zDim.name = @"z";
    GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
    GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
    
    GLFloat width = 1000;
    GLFloat height = 500;
    NSUInteger Nx = 4;
    NSUInteger Ny = 4;
    GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
    xDim.name = @"x";
    GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
    yDim.name = @"y";
    
    GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
    [internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude];
    
    [self compareSolution: internalModes toAnalyticalSolutionWithN2: N2_0 latitude: latitude depth: H];
}

- (void)testInternalWaveFullModesSpectralOptimized {
	GLFloat N2_0 = 1.69e-4;
	GLFloat rho0 = 1025;
	GLFloat g = 9.81;
	GLFloat latitude = 33.0;
	GLFloat H = 300;
	
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 64 domainMin: -H length: H]; zDim.name = @"z";
	GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
	GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
	
	GLFloat width = 10e3;
	GLFloat height = 5e3;
	NSUInteger Nx = 32;
	NSUInteger Ny = 32;
	GLDimension *xDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Nx domainMin: -width/2 length: width];
	xDim.name = @"x";
	GLDimension *yDim = [[GLDimension alloc] initDimensionWithGrid: kGLPeriodicGrid nPoints: Ny domainMin: -height/2 length: height];
	yDim.name = @"y";
	
	GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
	[internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[xDim, yDim, zDim] forLatitude: latitude maximumModes: 31];
	
	[self compareSolution: internalModes toAnalyticalSolutionWithN2: N2_0 latitude: latitude depth: H];
}


- (void)testPerformanceSpectral {
	GLFloat N2_0 = 1.69e-4;
	GLFloat rho0 = 1025;
	GLFloat g = 9.81;
	GLFloat latitude = 33.0;
	GLFloat H = 300;
	GLFloat k = 0.1;
	
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 256 domainMin: -H length: H]; zDim.name = @"z";
	GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
	GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
	
	GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];

	[self measureBlock:^{
        [internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: k forLatitude: latitude];
    }];
}

- (void)testPerformanceOptimized {
	GLFloat N2_0 = 1.69e-4;
	GLFloat rho0 = 1025;
	GLFloat g = 9.81;
	GLFloat latitude = 33.0;
	GLFloat H = 300;
	GLFloat k = 0.1;
	
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 256 domainMin: -H length: H]; zDim.name = @"z";
	GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
	GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
	
	GLInternalModesSpectral *internalModes = [[GLInternalModesSpectral alloc] init];
	[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: k forLatitude: latitude maximumModes: 0 zOutDim: zDim];
	
	[self measureBlock:^{
		[internalModes internalWaveModesFromDensityProfile: rho_bar wavenumber: k forLatitude: latitude maximumModes: 0 zOutDim: zDim];
	}];
}

@end
