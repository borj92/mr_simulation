function [ tfh ] = ReadTrajHeader(trajID)
%ReadTrajHeader reads the information stored in the 3-line header of a .traj file.
% tfh is a struct containing the duration, number of walkers and tmax. 
%
%
% Author: Ben Jordan (rmapbjo@ucl.ac.uk)

tfh.duration = fread(trajID,1,'double');
tfh.numwalkers = fread(trajID,1,'double');
tfh.tmax = fread(trajID,1,'double');


end

