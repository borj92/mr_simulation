% A script to run a batch of Camino simulations with different SNR.

%Move to correct working directory /Users/borj/Documents/PHD/MR Simulation Project

cd('/Users/borj/Documents/PHD/MR Simulation Project');
for cellrad = 5:12
    radius = cellrad * 1e-6;
    separation = num2str((2*radius) + 1e-6);
    for snr = 100:-1:1
        system(['/Users/borj/camino/bin/datasynth -walkers 100000 -tmax 1000 -snr ', num2str(snr), ' -voxels 1 -p 0.0 -schemefile v.scheme -initial uniform -substrate ply -plyfile sphere_rad',...
            num2str(cellrad), '.ply -meshsep ', separation, ' ', separation, ' ', separation, ' > ./Bfloat/plyfileout', num2str(cellrad), 'um_snr', num2str(snr), '.bfloat'])
    end
end




