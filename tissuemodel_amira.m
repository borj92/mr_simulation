function [tissue] = tissuemodel_amira(ficvf, fee, fvasc, radius, di, dv, amira)

dimx = max(amira(4).Val(:,1)) * 1e-6;
dimy = max(amira(4).Val(:,2)) * 1e-6;
dimz = max(amira(4).Val(:,3)) * 1e-6;

%Calculate voxel volume
voxvol = dimx*dimy*dimz;

%Calculate cell volume
cellvol = 4/3 * pi * (radius(1)^3);

%Use ficvf = 0.2 to ensure cell packing completes
ficvf2 = 0.2;

%Calculate number of cells to pack
numcells = floor((ficvf2*voxvol)/cellvol);

%Pack with cells
[centres, rads] = sampleSpheres([dimx, dimy, dimz], numcells, radius);
T = delaunayn(centres);

%Create Output

tissue.dimx = dimx;
tissue.dimy = dimy;
tissue.dimz = dimz;
tissue.ficvf = ficvf;
tissue.fee = fee;
tissue.fvasc = fvasc;
tissue.rads = rads;
tissue.voxvol = voxvol;
tissue.cellvol = cellvol;
tissue.numcells = numcells;
tissue.di = di;
tissue.dv = dv;
tissue.centres = centres;
tissue.T = T;
tissue.cellrad = radius;
tissue.adata = amira;

tissue;

end