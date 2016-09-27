%% SimToTraj is a matlab script to produce .traj files from Monte Carlo diffusion Simulation

cellrad = 8e-6;
tissue = tissuemodel(0.00047,0.00047,0.00047,0.3,0.7,[cellrad, 0], 1e-9);

for i = 1:10
    [traj, headerinfo] = tissuediffusion(tissue,1000,0.2,1e-4);
    WriteTraj(headerinfo, traj, ['./traj_files/mytraj_1000walkers_8umcellsonly_', num2str(i),'.traj']);
    matlabmail('rmapbjo@ucl.ac.uk',['mytraj_1000walkers_8umcellsonly_', num2str(i), '.traj'],'Simulation Complete', 'borjordanmatlab@gmail.com', 'matlabmail');
end