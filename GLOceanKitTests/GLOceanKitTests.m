//
//  GLOceanKitTests.m
//  GLOceanKitTests
//
//  Created by Jeffrey J. Early on 10/2/14.
//  Copyright (c) 2014 Jeffrey J. Early. All rights reserved.
//

#import <Cocoa/Cocoa.h>
#import <XCTest/XCTest.h>
#import <GLNumericalModelingKit/GLNumericalModelingKit.h>
#import "GLOceanKit.h"

@interface GLOceanKitTests : XCTestCase

@end

@implementation GLOceanKitTests

- (void)setUp {
    [super setUp];
    // Put setup code here. This method is called before the invocation of each test method in the class.
}

- (void)tearDown {
    // Put teardown code here. This method is called after the invocation of each test method in the class.
    [super tearDown];
}

- (void)testExample {
	GLFloat N2_0 = 1.69e-4;
	GLFloat rho0 = 1025;
	GLFloat g = 9.81;
	
	GLEquation *equation = [[GLEquation alloc] init];
	GLDimension *zDim = [[GLDimension alloc] initDimensionWithGrid: kGLEndpointGrid nPoints: 128 domainMin: -300.0 length: 300.0];
	GLFunction *z = [GLFunction functionOfRealTypeFromDimension:zDim withDimensions:@[zDim] forEquation:equation];
	GLFunction *rho_bar = [[z times: @(-N2_0*rho0/g)] plus: @(rho0)];
	
	GLInternalModes *internalModes = [[GLInternalModes alloc] init];
	[internalModes internalWaveModesFromDensityProfile: rho_bar withFullDimensions:@[zDim] forLatitude: 33.0];

    XCTAssert(YES, @"Pass");
}

- (void)testPerformanceExample {
    // This is an example of a performance test case.
    [self measureBlock:^{
        // Put the code you want to measure the time of here.
    }];
}

@end
