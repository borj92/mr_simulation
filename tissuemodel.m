%% Tissue Model sets up a single voxel and populates it with simulated cells.

function [tissue] = tissuemodel(dimx, dimy, dimz, ficvf, fee, radius, di)

%Check Inputs
if ficvf+fee ~= 1
    error('FICVF + FEE must == 1')
end

%Calculate Voxel volume
voxvol = dimx*dimy*dimz;

%Calculate cell volume
cellvol = 4/3 * pi * (radius(1)^3);

%Calculate number of cells
numcells = floor((ficvf*voxvol)/cellvol);

%Populate with spheres
[centres,rads] = sampleSpheres([dimx,dimy,dimz],numcells,radius);
T = delaunayn(centres);

%Create Output

tissue.dimx = dimx;
tissue.dimy = dimy;
tissue.dimz = dimz;
tissue.ficvf = ficvf;
tissue.fee =  fee;
tissue.rads = rads;
tissue.voxvol = voxvol;
tissue.cellvol = cellvol;
tissue.numcells = numcells;
tissue.di = di;
tissue.centres = centres;
tissue.T = T;
tissue.cellrad = radius;

tissue

end




