function initWithHorizontalWaveNumberSpectrum(self,GMAmplitude,options)
% initialize with a Alternative Interal Wave Spectrum in 
% function of horizontal wave number and mode
% 
% This only initializes the wave components, A0 is left untouched.
%
% - Topic: Initial conditions — Waves
% - Declaration: initWithHorizontalWaveNUmberSpectrum(GMAmplitude,options)
% - Parameter GMAmplitude: 
% - Parameter j_star: (optional) 
% - Parameter slope: (optional)


arguments
    self WVTransform {mustBeNonempty}
    GMAmplitude (1,1) double
    options.j_star (1,1) double = 3
    options.slope (1,1) double = 1
    options.shouldRandomizeAmplitude = 1    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a reasonable total wavenumber (Radial) axis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kh= self.Kh;
kRadial = self.radialWavenumberAxis();
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Distribution of Energy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j_star=options.j_star;
slope=options.slope;

% GM Parameters
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E_T = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;
%  E = E*(self.Lz/L_gm); % This correction fixes the amplitude so that the HKE variance at a given depth matches (instead of depth integrated energy)

% Compute the proper vertical function normalization
M = (j_star^2 +(2:1024).^2).^((-5/4));
M_norm = sum(M);

%Create the energy matrix 3D 
TotalEnergy = zeros(size(Kh));
Energy_slice=zeros(size(Kh(:,:,1)));

%%%% Redistributing the energy %%%
for j=(2:length(self.j))
    h=self.h(1,1,j);
    Kh_2D = Kh(:,:,j);
    LR= sqrt(self.g*h)/self.f;
    
    fun = @(k) (1./(k.^2*LR^2 + 1).^(1*slope))*LR;
    B_norm = integral(fun,kRadial(1),kRadial(end));

    for i=(1:length(kRadial)-1)       
         
        % Integrate the energy btw 2 Kh
        E = E_T*(integral(fun,kRadial(i),kRadial(i+1))/B_norm)*(((j^2 + j_star^2).^((-5/4)))/M_norm);
    
        % find all the kl point btw the two values of Kh
        ind = find(Kh_2D>=kRadial(i) & Kh_2D<kRadial(i+1));
    
        % Distribuite equally the energy btw all the points that are 
        % btw the circles (values of Kh)
        n = length(ind);
        if n > 0
            Energy_slice(ind) = E/n;
        end        
    end
   TotalEnergy(:,:,j)= Energy_slice;
    
    %clear the slice but maintain the size
    Energy_slice=zeros(size(Kh(:,:,1)));
end


disp([' Initial Total energy:',num2str(E_T), ', Total energy after distribution:',num2str(sum(TotalEnergy(:)))])
fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*sum(TotalEnergy(:))/E_T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% After comput the amplitude I insert that in the model 
% to get the variables in real space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The amplitude is:

shouldRandomizeAmplitude=options.shouldRandomizeAmplitude;

A = sqrt((TotalEnergy./self.h)/2);

 if shouldRandomizeAmplitude == 1 
        A_plus = A.*self.generateHermitianRandomMatrix( );
        A_minus = A.*self.generateHermitianRandomMatrix(  );

        %self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext).*randn( size(self.offgridModes.h_ext) );
        %self.offgridModes.PrecomputeExternalWaveCoefficients();                
    else
        % Randomize phases, but keep unit length
        A_plus = self.generateHermitianRandomMatrix(shouldExcludeNyquist=1, allowMeanPhase=1 );
        A_minus = self.generateHermitianRandomMatrix(shouldExcludeNyquist=1, allowMeanPhase=1 );

        goodIndices = abs(A_plus) > 0;
        A_plus(goodIndices) = A_plus(goodIndices)./abs(A_plus(goodIndices));
        A_plus = A.*A_plus;

        goodIndices = abs(A_minus) > 0;
        A_minus(goodIndices) = A_minus(goodIndices)./abs(A_minus(goodIndices));
        A_minus = A.*A_minus;

        % Check this factor of 2!!! Is the correct? squared
        % velocity to energy, I think.
        %self.offgridModes.U_ext = sqrt(2*GM3Dext./self.offgridModes.h_ext);
        %self.offgridModes.PrecomputeExternalWaveCoefficients();   
 end

self.Ap=A_plus;
self.Am=A_minus;

end