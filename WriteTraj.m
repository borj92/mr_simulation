function [arun_traj] = WriteTraj(headerinfo,traj,filename)
%WriteTraj writes the walker trajectories stored in the variable traj to the filename specified,
%after reordering the values in traj according to the specifications of the Camino diffusion toolkit
%software.
% 
% Inputs:
% 
% - headerinfo -> A structure containing the header information for the traj file. The three subfields are
%     - headerinfo.duration -> the duration of the simulation dynamics (0.2 for verdict)
%     - headerinfo.nwalkers -> the number of walkers in the simulation
%     - headerinfo.tmax     -> the maximum number of timesteps

%% Re-order traj matrix to conform to Camino software (time->index rather than index->time)

tmp = permute(traj,[3 1 2]);
arun_traj = reshape(tmp, (size(traj,1)*size(traj,3)), 6); %arun_traj is named in honour of Arun Niranjan -
% programming mastermind and general geezer.
%arun_traj has form -> | time | index | X | Y | Z | S | and is sorted in time-ascending order.

%% Create file and write header

%Create file using fopen.
fid = fopen(filename, 'w');

%Reopen file using specified settings

%Write header
fwrite(fid,headerinfo.duration,'double',0,'ieee-be');
fwrite(fid,headerinfo.nwalkers,'double',0,'ieee-be');
fwrite(fid,headerinfo.tmax,'double',0,'ieee-be');

%% Begin filling traj file with trajectories

for i = 1:size(arun_traj,1)
    
    %Write time (double), index (int), x (double), y (double), z (double).
    fwrite(fid,arun_traj(i,1),'double','ieee-be');
    fwrite(fid,arun_traj(i,2),'int','ieee-be');
    fwrite(fid,arun_traj(i,3),'double','ieee-be');
    fwrite(fid,arun_traj(i,4),'double','ieee-be');
    fwrite(fid,arun_traj(i,5),'double','ieee-be');
    
end

fclose(fid);

    
    