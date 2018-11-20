function prob_params = load_problem_params(runDIR)
% prob_params = load_problem_params(runDIR)
%
% NSF Internal Wave / Vortical Mode Analysis
% Script for reading problem_params file used to set up Winters Boussinesq 3D model
%
% Loads all variables from problem_params file.
% Input: runDIR should give relative or absolute path to rundirectory where /input/problem_params lives
%
% Output: structure variable 'prob_params' containing all variables set in problem_params file:
%   runlabel, do_nonlinear, do_second_scalar, vertical_coriolis
%   p1, p2
%   nx,ny,nz
%   dt,t0,tf
%   Lx,Ly,Lz
%   scalar_kind (for s1, s2 scalars)
%   do_forcing, do_immersed
%   g, (f0,beta), y_pivot, rho_0	!! NOTE: y_pivot is computed as pivot_point*Ly, with default as in read_userdata.f90
%   nu, kappa(s1, s2)
%   high_order_operator
%   dgrad, scalar_scale(s1, s2), u0
%   nparticles, particle_write_inc
%
% Written by Caterina Massidda, 3/12/2017
% Last modified by Caterina Massidda 6/1/2017

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Textscan problem_params and extract the first 'column' of the file containing key variables

prob_params_fileID = fopen(sprintf('%s/input/problem_params', runDIR));
prob_params_data = textscan(prob_params_fileID,'%s %*[^\n]');
fclose(prob_params_fileID);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create blank structure variable
prob_params = struct;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read variables inside problem_params
prob_params.runlabel = cell2mat(prob_params_data{1}(1));			% if run has special label

prob_params.do_nonlinear = strcmp('.TRUE.',cell2mat(prob_params_data{1}(2)));	% logical, solve nonlinear eqns of motion?
prob_params.do_second_scalar = strcmp('.TRUE.',cell2mat(prob_params_data{1}(3))); 	% logical, evolution of 2nd scalar field?
prob_params.vertical_coriolis = strcmp('.TRUE.',cell2mat(prob_params_data{1}(4))); 	% logical, use full rotation vector (non-traditional Coriolis terms?)

prob_params.p1= str2num(cell2mat(prob_params_data{1}(5)));			% MPI decomposition parameter p1, splits (nx)  
prob_params.p2 = str2num(cell2mat(prob_params_data{1}(6)));			% MPI decomposition parameter p2, splits (ny,nz)  p1*p2=numprocs=np
prob_params.nx = str2num(cell2mat(prob_params_data{1}(7)));			% nx    split by (p1)
prob_params.ny = str2num(cell2mat(prob_params_data{1}(8)));			% ny    split by (p1,p2
prob_params.nz = str2num(cell2mat(prob_params_data{1}(9)));			% nz    split by (p2)   
prob_params.dt = str2num(cell2mat(prob_params_data{1}(10)));			% dt       (s)      
prob_params.t0 = str2num(cell2mat(prob_params_data{1}(11)));			% t_start  (s)    
prob_params.tf = str2num(cell2mat(prob_params_data{1}(12)));			% t_end    (s)    
prob_params.Lx = str2num(cell2mat(prob_params_data{1}(13)));			% Lx       (m)
prob_params.Ly = str2num(cell2mat(prob_params_data{1}(14)));			% Ly       (m)
prob_params.Lz = str2num(cell2mat(prob_params_data{1}(15)));			% Lz       (m)

prob_params.scalar_kind(1) = cell2mat(prob_params_data{1}(16));			% char, defn of scalar s1 (r = Boussinesq density, i.e. rho)
prob_params.scalar_kind(2) = cell2mat(prob_params_data{1}(17));			% char, defn of scalar s2 (p =  passive tracer)

prob_params.do_forcing = strcmp('.TRUE.',cell2mat(prob_params_data{1}(18)));	% logical, call user forcing routine?
prob_params.do_immersed = strcmp('.TRUE.',cell2mat(prob_params_data{1}(19)));% logical, Is an immersed boundary specified?
prob_params.g = str2num(cell2mat(prob_params_data{1}(20)));			% gravity (m/s^2)

prob_params.f(1) = str2num(cell2mat(prob_params_data{1}(21)));			% f0  Coriolis parameter  (rad/s) 

% Note, fortran code reads next entry as a complex number, then parses beta and pivot point - here simply parse via Matlab 'deal' call.
eval(['[prob_params.f(2),prob_params.y_pivot]=deal',cell2mat(prob_params_data{1}(22)),';'])% Coriolis beta (1/(m s)), and pivot point for beta plane (m)
if prob_params.y_pivot==0							% set alternate pivot, as in read_userdata.f90
  prob_params.y_pivot = prob_params.Ly/2;
else
  prob_params.y_pivot = prob_params.y_pivot*prob_params.Ly;			% as per read_userdata.f90
end

prob_params.rho_0 = str2num(cell2mat(prob_params_data{1}(23)));			% rho_0  constant reference density (kg/m^3)
prob_params.nu = str2num(cell2mat(prob_params_data{1}(24)));			% nu     viscosity                  (m^2/s)                

prob_params.kappa(1) = str2num(cell2mat(prob_params_data{1}(25)));		% kappa  diffusivity for s1         (m^2/s) 
prob_params.kappa(2) = str2num(cell2mat(prob_params_data{1}(26)));		% kappa  diffusivity for s2         (m^2/s)  
prob_params.high_order_operator=strcmp('.TRUE.',cell2mat(prob_params_data{1}(27)));	% logical, use high order diffusion operators?
prob_params.dgrad = str2num(cell2mat(prob_params_data{1}(28)));			% characteristic val of drho/dz  (kg/m^4) 

prob_params.scalar_scale(1) = str2num(cell2mat(prob_params_data{1}(29)));	% characteristic val of s1 fluctuations
prob_params.scalar_scale(2) = str2num(cell2mat(prob_params_data{1}(30)));	% characteristic val of s2 fluctuations 
prob_params.u0 = str2num(cell2mat(prob_params_data{1}(31)));			% characteristic velocity  (m/s)        
prob_params.nparticles = str2num(cell2mat(prob_params_data{1}(32)));		% number of particles satisfying Lagrangian eqn Dx/Dt 
prob_params.particle_write_inc = str2num(cell2mat(prob_params_data{1}(33)));	% number of time steps btwn writing particle data

