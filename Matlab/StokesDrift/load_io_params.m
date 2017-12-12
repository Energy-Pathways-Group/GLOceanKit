 function io_params=load_io_params(runDIR)
  % io_params = load_io_params(runDIR)
  %
  % NSF Internal Wave / Vortical Mode Analysis
  % Script for reading io_params file used to set up Winters Boussinesq 3D model
  %
  % Loads all variables from io_params file.
  % Input: runDIR should give relative or absolute path to rundirectory where /input/load_io_params lives
  % Output: structure variable 'io_params' containing all variables set in io_params file:
  %   num_file_sets, filename_root
  %   mode
  %   nsteps
  %   ilocs, jlocs, klocs
  %   variable_key, write_s1_bar, write_s2_bar 
  %
  % Written by Caterina Massidda, 5/12/2017
  % last modification 5/31/2017
  %
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  % Scan io_params and extract the first column of the file, containing key variables
  
  io_params_fileID = fopen(sprintf('%s/codes_etc/input/io_params', runDIR));
  io_params_data = textscan(io_params_fileID,'%s %*[^\n]');
  fclose(io_params_fileID);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %create blank structure variable
  io_params = struct;

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Read variables inside io_params
  io_params.num_file_sets = str2num(cell2mat(io_params_data{1}(1)));                  % num of filesets to process

  ind = 0;

  for j = 1:io_params.num_file_sets
    io_params.filename_root(j) = (io_params_data{1}(2+ind));                          % rootname of each file set (cell array)
    io_params.mode(j) = (io_params_data{1}(3+ind));                                   % append w/each write, or create new file
    io_params.nsteps(j) = str2num(cell2mat(io_params_data{1}(4+ind)));                % steps between writes  
    io_params.ilocs(j,1) = str2num(cell2mat(io_params_data{1}(5+ind)));               % i start
    io_params.ilocs(j,2) = str2num(cell2mat(io_params_data{1}(6+ind)));               % i end
    io_params.ilocs(j,3) = str2num(cell2mat(io_params_data{1}(7+ind)));               % i inc
    io_params.jlocs(j,1) = str2num(cell2mat(io_params_data{1}(8+ind)));               % j start
    io_params.jlocs(j,2) = str2num(cell2mat(io_params_data{1}(9+ind)));               % j end
    io_params.jlocs(j,3) = str2num(cell2mat(io_params_data{1}(10+ind)));              % j inc
    io_params.klocs(j,1) = str2num(cell2mat(io_params_data{1}(11+ind)));              % k start
    io_params.klocs(j,2) = str2num(cell2mat(io_params_data{1}(12+ind)));              % k end
    io_params.klocs(j,3) = str2num(cell2mat(io_params_data{1}(13+ind)));              % k inc
    io_params.variable_key(j,1) = str2num(cell2mat(io_params_data{1}(14+ind)));       % u variable_key is (5xio_params.num_file_sets) row1 u v w s1 and s2 fle1 
    io_params.variable_key(j,2) = str2num(cell2mat(io_params_data{1}(15+ind)));       % v
    io_params.variable_key(j,3) = str2num(cell2mat(io_params_data{1}(16+ind)));       % w
    io_params.variable_key(j,4) = str2num(cell2mat(io_params_data{1}(17+ind)));       % s1
    io_params.variable_key(j,5) = str2num(cell2mat(io_params_data{1}(18+ind)));       % s2
    io_params.write_s1_bar(j) = str2num(cell2mat(io_params_data{1}(17+ind)));         % write s1 (don't write if 0, write s1' if 1, s1' + s1_bar if 2)
    io_params.write_s2_bar(j) = str2num(cell2mat(io_params_data{1}(18+ind)));         % write s2 (don't write if 0, write s2' if 1, s2' + s2_bar if 2)

    ind = ind + 17;                                                                   % 17 is the number of parameters that can be set for each file
  end

end