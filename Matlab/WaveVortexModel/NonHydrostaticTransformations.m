        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            % hydrostatic modes commute with the DFT
%             u = permute(u,[3 1 2]); % keep adjacent in memory
%             u = reshape(u,self.Nz,[]);
%             u_bar = self.PF*u;
%             u_bar = reshape(u_bar,self.nModes,self.Nx,self.Ny);
%             u_bar = permute(u_bar,[2 3 1]);
%             u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2);

            u_bar = fft(fft(u,self.Nx,1),self.Ny,2);
            u_bar = permute(u_bar,[3 1 2]);
            u_bar = reshape(u_bar,self.Nz,[]);
            u_bar = self.PF*u_bar;
            u_bar = reshape(u_bar,self.nModes,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            % hydrostatic modes commute with the DFT
%             w = permute(w,[3 1 2]); % keep adjacent in memory
%             w = reshape(w,self.Nz,[]);
%             w_bar = self.QG*w;
%             w_bar = reshape(w_bar,self.nModes,self.Nx,self.Ny);
%             w_bar = permute(w_bar,[2 3 1]);
%             w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2);

            w_bar = fft(fft(w,self.Nx,1),self.Ny,2);
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nz,[]);
            w_bar = self.QG*w_bar;
            w_bar = reshape(w_bar,self.nModes,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            
        end
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
            % hydrostatic modes commute with the DFT
%             u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
%             u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
%             u_bar = reshape(u_bar,self.nModes,[]);
%             u = self.PFinv*u_bar;
%             u = reshape(u,self.Nz,self.Nx,self.Ny);
%             u = permute(u,[2 3 1]);

            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.nModes,[]);
            u_bar = self.PFinv*u_bar;
            u_bar = reshape(u_bar,self.Nz,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
        end
                
        function w = TransformToSpatialDomainWithG(self, w_bar )
            % hydrostatic modes commute with the DFT
%             w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
%             w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
%             w_bar = reshape(w_bar,self.nModes,[]);
%             w = self.QGinv*w_bar;
%             w = reshape(w,self.Nz,self.Nx,self.Ny);
%             w = permute(w,[2 3 1]);
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.nModes,[]);
            w_bar = self.QGinv*w_bar;
            w_bar = reshape(w_bar,self.Nz,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
        end
        
        function [u,ux,uy,uz] = TransformToSpatialDomainWithFAllDerivatives(self, u_bar)
%             u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
% 
%             u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
%             u_bar = reshape(u_bar,self.nModes,[]);
%             u = self.PFinv*u_bar;
%             u = reshape(u,self.Nz,self.Nx,self.Ny);
%             u = permute(u,[2 3 1]);
% 
%             ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
%             uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
% 
%             uz = self.QGinv*( squeeze(self.Q./self.P).*u_bar );
%             uz = reshape(uz,self.Nz,self.Nx,self.Ny);
%             uz = permute(uz,[2 3 1]);
%             uz = (-shiftdim(self.N2,-2)/self.g).*uz;
            
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.nModes,[]);
            u_bar = self.PFinv*u_bar;
            u_bar = reshape(u_bar,self.Nz,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);

            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');

            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
            
            u_bar = permute(u,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.Nz,[]);
            u_bar = self.PF*u_bar;
            uz = self.QGinv*( squeeze(self.Q./self.P).*u_bar );
            uz = reshape(uz,self.Nz,self.Nx,self.Ny);
            uz = permute(uz,[2 3 1]);
            uz = (-shiftdim(self.N2,-2)/self.g).*uz;
        end  
        
        function [w,wx,wy,wz] = TransformToSpatialDomainWithGAllDerivatives(self, w_bar )
%             w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
% 
%             w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
%             w_bar = reshape(w_bar,self.nModes,[]);
%             w = self.QGinv*w_bar;
%             w = reshape(w,self.Nz,self.Nx,self.Ny);
%             w = permute(w,[2 3 1]);
% 
%             wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
%             wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
%             
%             wz = self.PFinv* ( squeeze(self.P./(self.Q .* self.h)) .* w_bar);
%             wz = reshape(wz,self.Nz,self.Nx,self.Ny);
%             wz = permute(wz,[2 3 1]);

            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.nModes,[]);
            w_bar = self.QGinv*w_bar;
            w_bar = reshape(w_bar,self.Nz,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);

            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');

            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
            
            w_bar = permute(w,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.Nz,[]);
            w_bar = self.QG*w_bar;
            wz = self.PFinv* ( squeeze(self.P./(self.Q .* self.h)) .* w_bar);
            wz = reshape(wz,self.Nz,self.Nx,self.Ny);
            wz = permute(wz,[2 3 1]);
        end