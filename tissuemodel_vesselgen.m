%% Tissue Model sets up a single voxel and populates it with simulated cells and cylinders.

function [tissue] = tissuemodel_vesselgen(dimx, dimy, dimz, ficvf, fee, fvasc, radius, di, dv, vesseldensity)

%Check Inputs
% if ficvf+fee+fvasc > 1
%     error('FICVF + FEE + FVASC must == 1')
% end

%Calculate Voxel volume
voxvol = dimx*dimy*dimz;

%Calculate cell volume
cellvol = 4/3 * pi * (radius(1)^3);

%Use ficvf = 0.3 to calculate cell number
ficvf2 = 0.2; %BJ edit 14/6/16

%Calculate number of cells
numcells = floor((ficvf2*voxvol)/cellvol);

%Populate with spheres
[centres,rads] = sampleSpheres([dimx,dimy,dimz],numcells,radius);
T = delaunayn(centres);

%Populate with vessels
X = rand(vesseldensity(1),1) * dimx;
Y = rand(vesseldensity(2),1) * dimy;
Z = rand(vesseldensity(3),1) * dimz;
vtree = MST_tree(1,[0;X], [0;Y], [0;Z],[],[],[],[],'-b');
vtree = resample_tree(vtree, dimx/100);
vtree.D = vtree.D*2e-6;
child = child_tree(vtree);
route = PL_tree(vtree);
nodes = [vtree.X, vtree.Y, vtree.Z];
node_T = delaunayn(nodes);

%Create Output

tissue.dimx = dimx;
tissue.dimy = dimy;
tissue.dimz = dimz;
tissue.ficvf = ficvf;
tissue.fee =  fee;
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
tissue.vtree = vtree;
tissue.child = child;
tissue.vroute = route;
tissue.nodes = nodes;
tissue.node_T = node_T;

tissue

end




