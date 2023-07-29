ncfile = NetCDFFile('QGMonopole.nc');

%% Highly recommended ACDD conventions
% Descriptions of these four attributes are found here:
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Highly_Recommended
ncfile.addAttribute('title','CanonicalEddy'); 
ncfile.addAttribute('summary','A numerical model of a propagating, isolated, 1.5 layer eddy under quasigeostrophic balance on a beta plane, with trajectories of advected particles.'); 
%%Using this vocabulary for keywords as recommended by ACDD: https://gcmd.earthdata.nasa.gov/KeywordViewer/scheme/all?gtm_search=wave&gtm_scheme=all
ncfile.addAttribute('keywords','ocean currents, component process models, mesoscale eddies, potential vorticity, ocean waves, turbulence')
ncfile.addAttribute('Conventions','CF-1.10, ACDD-1.3'); 

%% Recommended ACDD conventions
% https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3#Recommended
ncfile.addAttribute('id','canonical_eddy'); %Could also be a DOI once we have that 
ncfile.addAttribute('naming_authority','com.jeffreyearly');%Could also be Zenodo 
ncfile.addAttribute('product_version','1.0.0');%suggested
ncfile.addAttribute('history',[datestr(now,1) 'v. 1.0.0 initial version']);
ncfile.addAttribute('source','WaveVortexModel version XXXX');
%This is the only CF convention that is not also an ACDD convention...seems redundant with creator_institution below but whatever 
ncfile.addAttribute('institution','NorthWest Research Associates');
ncfile.addAttribute('processing_level','Unprocessed original model output'); %textual description 
ncfile.addAttribute('comment','Initially created in support of ''A tensor-valued integral theorem for the gradient of a vector field, with a fluid dynamical application'', by Lilly, Feske, Fox-Kemper, and Early (2023).');
ncfile.addAttribute('acknowledgment','J. Early''s work on this model run was supported by grant number 2049521 from the Physical Oceanography program of the United States National Science Foundation.  Additional support for the model development was provided by XXX.');
ncfile.addAttribute('license','Creative Commons Attribution 4.0 International, https://creativecommons.org/licenses/by/4.0/');
ncfile.addAttribute('standard_name_vocabulary','CF Standard Name Table');

%in case you change in the future: 
%ncfile.addAttribute('date_modified',datestr(now,1));%suggested
ncfile.addAttribute('creator_name','Jeffrey J. Lilly');
ncfile.addAttribute('creator_email','jearly@nwra.com');
ncfile.addAttribute('creator_url','https://jeffreyearly.com');
ncfile.addAttribute('creator_institution','NorthWest Research Associates'); %suggested
ncfile.addAttribute('project','Global eddy-driven transport estimated from in situ Lagrangian observations');
ncfile.addAttribute('publisher_name','Jonathan M. Lilly') %The name of the person (or other entity specified by the publisher_type attribute) responsible for publishing the data file or product to users, with its current metadata and format.
ncfile.addAttribute('publisher_email','jmlilly@psi.edu');
ncfile.addAttribute('publisher_url','http://www.jmlilly.net')
ncfile.addAttribute('publisher_institution','Planetary Science Institute') %suggested
%skipping these unless you want to acknowledge another contributor 
%ncfile.addAttribute('contributor_name',''); %suggested
%ncfile.addAttribute('contributor_role',''); %suggested
%skipping these recommended attributes as they seem like overkill:
%geospatial_bounds, geospatial_bounds_crs, geospatial_bounds_vertical_crs
ncfile.addAttribute('geoespatial_lat_min',min(lat));%XXX
ncfile.addAttribute('geoespatial_lat_max',max(lat));%XXX
ncfile.addAttribute('geoespatial_lat_units','degree_north');
ncfile.addAttribute('geoespatial_lon_min',min(lon));%XXX
ncfile.addAttribute('geoespatial_lon_max',max(lon));%XXX
ncfile.addAttribute('geoespatial_lon_units','degree_east');
ncfile.addAttribute('time_coverage_start',datestr(min(num)));%XXX
ncfile.addAttribute('time_coverage_end',datestr(max(num)));%XXX
ncfile.addAttribute('time_coverage_resolution','1h');%XXX
ncfile.addAttribute('platform','Models'); %suggested
ncfile.addAttribute('platform_vocabulary','GCMD, https://gcmd.earthdata.nasa.gov/static/kms/'); %suggested
%From https://docs.unidata.ucar.edu/netcdf-java/4.6/userguide/metadata/DataDiscoveryAttConvention.html
%The "cdm_data_type" attribute gives the THREDDS data type appropriate for this dataset. E.g., "Grid", "Image", "Station", "Trajectory", "Radial".
ncfile.addAttribute('cdm_data_type','Grid'); 
ncfile.addAttribute('references','Early, Lelong, and Sundermeyer (2022).  A generalized wave-vortex decomposition for rotating Boussinesq flows with arbitrary stratification, https://doi.org/10.1017/jfm.2020.995. Lilly, Feske, Fox-Kemper, and Early (2023).  A tensor-valued integral theorem for the gradient of a vector field, with a fluid dynamical application, submitted.'


ncfile.close();