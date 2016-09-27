function [tfl] = ReadTrajLine(trajID)
%A function to read a line from a trajfile and return the walker info.
%
%
% Author: Ben Jordan (rmapbjo@ucl.ac.uk)

tfl.time = fread(trajID,1,'double');
tfl.index = fread(trajID,1,'int');
tfl.x = fread(trajID,1,'double');
tfl.y = fread(trajID,1,'double');
tfl.z = fread(trajID,1,'double');



end

