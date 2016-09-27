function [fid] = OpenTraj(filename)
%OpenTraj opens a traj file using the painstakingly obtained machine format and encoding settings
%
%
% Author: Ben Jordan (rmapbjo@ucl.ac.uk)

fid = fopen(filename,'rb','ieee-be','ISO-8859-1');


end

