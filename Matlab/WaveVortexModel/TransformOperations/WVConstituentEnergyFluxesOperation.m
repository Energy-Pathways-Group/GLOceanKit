classdef WVConstituentEnergyFluxesOperation < WVOperation
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Property1
    end

    methods
        function self = WVConstituentEnergyFluxesOperation()
            outputVar = WVVariableAnnotation.empty;
            outputVar(1) = WVVariableAnnotation('E0_all_all',{},'m^3/s^2', 'flux into all, from all');
            outputVar(end+1) = WVVariableAnnotation('E0_igw_igw',{},'m^3/s^2', 'flux into E0, from wave grad-wave');
            outputVar(end+1) = WVVariableAnnotation('E0_igw_io',{},'m^3/s^2', 'flux into E0, from wave grad-inertial');
            outputVar(end+1) = WVVariableAnnotation('E0_igw_g',{},'m^3/s^2', 'flux into E0, from wave grad-geostrophic');
            outputVar(end+1) = WVVariableAnnotation('E0_io_igw',{},'m^3/s^2', 'flux into E0, from inertial grad-wave');
            outputVar(end+1) = WVVariableAnnotation('E0_io_g',{},'m^3/s^2', 'flux into E0, from inertial grad-geostrophic');
            outputVar(end+1) = WVVariableAnnotation('E0_g_igw',{},'m^3/s^2', 'flux into E0, from geostrophic grad-wave');
        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end