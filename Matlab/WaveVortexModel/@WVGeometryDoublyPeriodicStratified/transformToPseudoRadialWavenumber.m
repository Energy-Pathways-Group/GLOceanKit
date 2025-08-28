function [varargout] = transformToPseudoRadialWavenumber(self,energyReservoir,varargin)     
% transforms in the from (j,kRadial) to kPseudoRadial
%
% Sums all the variance/energy in radial bins `kPseudoRadial`.
%
% - Topic: Operations — Transformations
% - Declaration: [varargout] = transformToRadialWavenumber(varargin) 
% - Parameter varargin: variables with dimensions $$(j,kl)$$
% - Returns varargout: variables with dimensions $$(kRadial)$$ or $$(kRadial,j)$$

% Thi is the final output axis for wavenumber
varargout = cell(size(varargin));
switch energyReservoir
    case EnergyReservoir.geostrophic_kinetic
    case EnergyReservoir.geostrophic_kinetic_mda
    case EnergyReservoir.geostrophic_potential
    case EnergyReservoir.geostrophic_potential_mda
    case EnergyReservoir.geostrophic
    case EnergyReservoir.mda
    case EnergyReservoir.geostrophic_mda
        varargout{:} = self.transformToPseudoRadialWavenumberA0(varargin{:});
    case EnergyReservoir.igw
    case EnergyReservoir.io
    case EnergyReservoir.wave
        varargout{:} = self.transformToPseudoRadialWavenumberApm(varargin{:});
    case EnergyReservoir.total
        if self.wvt.isHydrostatic
            varargout{:} = self.transformToPseudoRadialWavenumberA0(varargin{:});
        else
            error("The total energy cannot be combined in this way for a non-hydrostatic simulation.");
        end
    otherwise
        error("unknown energy reservoir");
end

end