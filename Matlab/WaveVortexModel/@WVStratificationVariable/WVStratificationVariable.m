classdef WVStratificationVariable < WVStratification & CAAnnotatedClass
    properties (Access=public) %(GetAccess=private, SetAccess=private) %(Access=private)
        dLnN2
        
        % Transformation matrices
        PF0inv, QG0inv % size(PFinv,PGinv)=[Nz x Nj]
        PF0, QG0 % size(PF,PG)=[Nj x Nz]
        h_0 % [Nj 1]

        P0 % Preconditioner for F, size(P)=[Nj 1]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat
        Q0 % Preconditioner for G, size(Q)=[Nj 1]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat.

        zInterp
        PFinvInterp, QGinvInterp
    end

    properties (Dependent)
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
        Lr2
    end

    methods
        function self = WVStratificationVariable(Lz, Nz, options, directInit)
            % create matrices for hydrostatics in variable stratification
            %
            % To initialize:
            % 1) Pass (Lz,Nz,N2) as a minimum set, and defaults will be
            % chosen (or rho instead of N2).
            % 2) Pass (Lz,Nz,N2,Nj,latitude,rho0) to fully, uniquely
            % specify the transforms
            % 3) 
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter rho:  (optional) function_handle specifying the density as a function of depth on the domain [-Lz 0]
            % - Parameter stratification:  (optional) function_handle specifying the stratification as a function of depth on the domain [-Lz 0]
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Parameter rho0: (optional) density at the surface z=0 (default is 1025 kg/m^3)
            % - Returns wvt: a new WVTransformHydrostatic instance
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.j (:,1) double {mustBeNonempty}
                options.Nj (1,1) double {mustBePositive}
                options.rhoFunction function_handle = @isempty
                options.N2Function function_handle = @isempty
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81

                % ALL of these must be set for direct initialization to
                % avoid actually computing the modes.
                directInit.dLnN2 (:,1) double
                directInit.PF0inv
                directInit.QG0inv
                directInit.PF0
                directInit.QG0
                directInit.h_0 (:,1) double
                directInit.P0 (:,1) double
                directInit.Q0 (:,1) double
                directInit.z_int (:,1) double
            end
            superclassOptions = namedargs2cell(options);
            self@WVStratification(Lz,Nz,superclassOptions{:});
            allFields = cell2struct([struct2cell(options);struct2cell(directInit)],[fieldnames(options);fieldnames(directInit)]);
            canInitializeDirectly = all(isfield(allFields, WVStratificationVariable.namesOfRequiredPropertiesForStratification));

            if canInitializeDirectly == true
                self.dLnN2 = directInit.dLnN2;
                self.PF0inv = directInit.PF0inv;
                self.QG0inv = directInit.QG0inv;
                self.PF0 = directInit.PF0;
                self.QG0 = directInit.QG0;
                self.P0 = directInit.P0;
                self.Q0 = directInit.Q0;
                self.h_0 = directInit.h_0;
                self.z_int = directInit.z_int;
            else
                self.dLnN2 = self.verticalModes.rho_zz./self.verticalModes.rho_z;
                [self.P0,self.Q0,self.PF0inv,self.PF0,self.QG0inv,self.QG0,self.h_0,self.z_int] = self.verticalProjectionOperatorsForGeostrophicModes(self.Nj);
            end

            % self.dftBuffer = zeros(self.spatialMatrixSize);
            % self.wvBuffer = zeros([self.Nz self.Nkl]);
            % [self.dftPrimaryIndex, self.dftConjugateIndex, self.wvConjugateIndex] = self.horizontalModes.indicesFromWVGridToDFTGrid(self.Nz,isHalfComplex=1);
        end

        % function initializeStratifiedFlow(wvt)
        %     % After initializing the WVTransform, this method can be called
        %     % and the WVStratifiedFlow will register.
        %     arguments
        %         wvt WVTransform
        %     end
        %     wvt.addPropertyAnnotations(WVStratificationVariable.propertyAnnotationsForStratifiedFlow);
        %     wvt.addOperation(EtaTrueOperation());
        %     wvt.addOperation(APVOperation());
        % end

        function Lr2 = get.Lr2(self)
            Lr2 = self.g*self.h_0/(self.f*self.f);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformation matrices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Finv = get.FinvMatrix(wvt)
            % transformation matrix $$F^{-1}$$
            %
            % A matrix that transforms a vector from vertical mode space to physical
            % space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Finv = FinvMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVStratification
            end
            Finv = shiftdim(wvt.P0,1) .* wvt.PF0inv;
        end

        function F = get.FMatrix(wvt)
            % transformation matrix $$F$$
            %
            % A matrix that transforms a vector from physical
            % space to vertical mode space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: F = FMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVStratification
            end
            F = wvt.PF0 ./ shiftdim(wvt.P0,2);
        end

        function Ginv = get.GinvMatrix(wvt)
            % transformation matrix $$G^{-1}$$
            %
            % A matrix that transforms a vector from vertical mode space to physical
            % space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Ginv = GinvMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVStratification
            end
            Ginv = shiftdim(wvt.Q0,1) .* wvt.QG0inv;
        end

        function G = get.GMatrix(wvt)
            % transformation matrix $$G$$
            %
            % A matrix that transforms a vector from physical
            % space to vertical mode space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: G = GMatrix(wvt)
            % - Returns Ginv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVStratification
            end
            G = wvt.QG0 ./ shiftdim(wvt.Q0,2);
        end
    end

    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVStratificationVariable.propertyAnnotationsForStratification();
        end

        function vars = classRequiredPropertyNames()
            vars = WVStratificationVariable.namesOfRequiredPropertiesForStratification();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForStratification()
            requiredPropertyNames = WVStratification.namesOfRequiredPropertiesForStratification();
            requiredPropertyNames = union(requiredPropertyNames,WVStratificationVariable.newRequiredPropertyNames());
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'dLnN2','PF0inv','QG0inv','PF0','QG0','P0','Q0','h_0','z_int'};
        end

        function propertyAnnotations = propertyAnnotationsForStratification()
            % return array of WVPropertyAnnotation initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratificationVariable class (as well as its
            % superclass).
            %
            % - Topic: Internal
            % - Declaration: propertyAnnotations = WVStratificationVariable.propertyAnnotationsForStratifiedFlow()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = WVStratification.propertyAnnotationsForStratification();
            propertyAnnotations(end+1) = CANumericProperty('PF0inv',{'z','j'},'','Preconditioned F-mode inverse transformation');
            propertyAnnotations(end+1) = CANumericProperty('QG0inv',{'z','j'},'','Preconditioned G-mode inverse transformation');
            propertyAnnotations(end+1) = CANumericProperty('PF0',{'j','z'},'','Preconditioned F-mode forward transformation');
            propertyAnnotations(end+1) = CANumericProperty('QG0',{'j','z'},'','Preconditioned G-mode forward transformation');
            propertyAnnotations(end+1) = CANumericProperty('P0',{'j'},'','Preconditioner for F, size(P)=[1 Nj]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat');
            propertyAnnotations(end+1) = CANumericProperty('Q0',{'j'},'','Preconditioner for G, size(Q)=[1 Nj]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. ');

            propertyAnnotations(end+1) = CANumericProperty('h_0',{'j'},'m', 'equivalent depth of each geostrophic mode', detailedDescription='- topic: Domain Attributes — Stratification');
            propertyAnnotations(end+1) = CANumericProperty('Lr2',{'j'},'m^2', 'squared Rossby radius');
        end

        function [Lz,Nz,options] = requiredPropertiesForStratificationFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options
            end
            [Lz,Nz,stratOptions] = WVStratification.requiredPropertiesForStratificationFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVStratificationVariable.newRequiredPropertyNames);
            newOptions = namedargs2cell(vars);
            options = cat(2,stratOptions,newOptions);
        end

        function stratification = stratificationFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                stratification WVStratificationVariable {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            stratification = WVStratificationVariable.stratificationFromGroup(ncfile);
        end

        function stratification = stratificationFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                stratification WVStratificationVariable {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVStratificationVariable.namesOfRequiredPropertiesForStratification);
            [Lz,Nz,options] = WVStratificationVariable.requiredPropertiesForStratificationFromGroup(group);
            stratification = WVStratificationVariable(Lz,Nz,options{:});
        end
        
    end
end