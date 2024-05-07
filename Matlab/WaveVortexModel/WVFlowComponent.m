classdef WVFlowComponent < handle
    %Orthogonal solution group
    %
    % Each degree-of-freedom in the model is associated with an analytical
    % solution to the equations of motion. This class groups together
    % solutions of a particular type and provides a mapping between their
    % analytical solutions and their numerical representation.
    %
    % Perhaps the most complicate part of the numerical implementation is
    % the indexing---finding where each solution is represented
    % numerically. In general, a solution will have some properties, e.g.,
    %   (kMode,lMode,jMode,phi,A,omegasign) 
    % which will have a primary and conjugate part, each of which might be
    % in two different matrices.
    %
    % - Topic: Initialization
    properties (Access=private)
        bitmask = 0
    end
    properties
        % name of the flow feature
        %
        % long-form version of the feature name, e.g., "internal gravity wave"
        % - Topic: Properties
        name

        % name of the flow feature
        %
        % camel-case version of the feature name, e.g., "internalGravityWave"
        % - Topic: Properties
        camelCaseName

        % abbreviated name
        %
        % abreviated feature name, e.g., "igw" for internal gravity waves.
        % - Topic: Properties
        abbreviatedName

        % reference to the wave vortex transform
        %
        % reference to the WVTransform instance
        % - Topic: Properties
        wvt

        % returns a mask indicating where solutions live in the Ap matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % different solution types live in the Ap matrix.
        %
        % - Topic: Masks
        maskAp

        % returns a mask indicating where solutions live in the Am matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % different solution types live in the Am matrix.
        %
        % - Topic: Masks
        maskAm

        % returns a mask indicating where solutions live in the A0 matrix.
        %
        % Returns a 'mask' (matrix with 1s or 0s) indicating where
        % different solution types live in the A0 matrix.
        %
        % - Topic: Masks
        maskA0
    end
    methods
        function self = WVFlowComponent(wvt,options)
            % create a new orthogonal solution group
            %
            % - Topic: Initialization
            % - Declaration:  solnGroup = WVFlowComponent(wvt)
            % - Parameter wvt: instance of a WVTransform
            % - Returns solnGroup: a new orthogonal solution group instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.maskAp = 0
                options.maskAm = 0
                options.maskA0 = 0
            end
            self.wvt = wvt;
            self.maskAp = options.maskAp;
            self.maskAm = options.maskAm;
            self.maskA0 = options.maskA0;
        end

    end
    
end

