function [traj,header] = ReadTraj(filename,numlines)
%ReadTraj is a function to read the contents of a trajectory file, in order to check they have been
%written correctly.
%
%
% Author: Ben Jordan (rmapbjo@ucl.ac.uk)

fid = OpenTraj(filename);

header = ReadTrajHeader(fid);

for i = 1:numlines
    
    traj(i) = ReadTrajLine(fid);
    
end

fclose(fid);

end