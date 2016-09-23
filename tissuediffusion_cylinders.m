function [traj, headerinfo, stats] = tissuediffusion_cylinders(tissue, numwalkers, max_t, dt, initialisation)
% Monte Carlo simulation of diffusion MR acquisition from water molecule diffusion through tissue.
%
% Input: tissue - A structure containing tissue properties:
%  
%             tissue.dimx     -   Voxel X-dimension
%             tissue.dimy     -   Voxel Y-dimension
%             tissue.dimz     -   Voxel Z-dimension
%             tissue.ficvf    -   Intracellular Volume Fraction
%             tissue.fee      -   Extracellular Volume Fraction
%             tissue.rads     -   Column-vector of cell radii
%             tissue.voxvol   -   Total Volume of Voxel
%             tissue.cellvol  -   Total Volume of all cells
%             tissue.numcells -   Number of cells in voxel
%             tissue.labels   -   Label matrix of tissue types
%             tissue.centres  -   N*3 matrix of cell positions
%             tissue.di       -   Diffusivity of extracellular space
%          
%         scheme - A matrix containing pulse sequence gradients
%
%         numwalkers - Number of walkers to be propagated
%
%         max_t - maximum simulation time per walker (s)
%
%         dt - Time per step (s)
% 
%         initialisation - String specifying how walkers are to be initialised.
%             -> 'random_full'        - Random distribution across whole voxel.
%             -> 'random_quartile'    - Random distribution across inner half. DEFAULT
%             -> 'random_cells'       - Only initialise in cells.
%             -> 'fixed_cells'        - Initialise within cells in the order they are saved.
%             -> 'vasc'               - Only initialise in vasculature. USE SMALL NUMWALKERS, TAKES
%             BLOOMIN AGES.
%             -> 'free_full'          - Only initialise in EES, full voxel.
%             -> 'free_quartile       - Only initialise in EES, inner half.
% 
% Outputs: traj - Matrix of walker trajectories
%          headerinfo - Structure containing information needed to create trajectory files
%
% Author: Ben Jordan - rmapbjo@ucl.ac.uk

%-----------------------------------------------------------------------------
% Funtions for walker propagation
%-----------------------------------------------------------------------------

%Initialise walker properties
    function [walker,stats] = setwalker(stats, initialisation)
        
    % Initialise walker positions based on user input 'initialisation'
    switch initialisation
        case 'random_full'
            walker.x = rand * tissue.dimx;
            walker.y = rand * tissue.dimy;
            walker.z = rand * tissue.dimz;
            
        case 'random_quartile'
            chance = rand;
            if(chance < tissue.ficvf)
                index = ceil(rand * length(tissue.centres));
                walker.x = tissue.centres(index,1);
                walker.y = tissue.centres(index,2);
                walker.z = tissue.centres(index,3);
                walker.lab = 1;
                
            elseif(chance >= tissue.ficvf && chance < (tissue.ficvf + tissue.fee))                
                walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
                walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
                walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
                walker.lab = NaN;
                walker = isintra(walker,tissue);
                while(walker.lab~=0)
                    walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
                    walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
                    walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
                    walker = isintra(walker,tissue);
                end
                
            else
                walker.x = tissue.vtree.X(2);
                walker.y = tissue.vtree.Y(2);
                walker.z = tissue.vtree.Z(2);
                walker.lab = 2;
                
            end
            
        case 'random_cells'
            %Initialise only in intercellular space, random cell choice
            index = ceil(rand * length(tissue.centres));
            walker.x = tissue.centres(index,1);
            walker.y = tissue.centres(index,2);
            walker.z = tissue.centres(index,3);
            
        case 'fixed_cells'
            %Initialise only in intercellular space, systematic cell choice
            walker.x = tissue.centres(i,1);
            walker.y = tissue.centres(i,2);
            walker.z = tissue.centres(i,3);
            
        case 'vasc'
            %Initialise only in vasculature
            walker.x = tissue.vtree.X(1);
            walker.y = tissue.vtree.Y(1);
            walker.z = tissue.vtree.Z(1);
            walker.lab = 2;
            
        case 'random_quartile_novasc'
            if(tissue.fvasc ~= 0)
                error('Specified initialisation requires fvasc = 0');
            end
            chance = rand;
            if(chance < tissue.ficvf)
                index = ceil(rand * length(tissue.centres));
                walker.x = tissue.centres(index,1);
                walker.y = tissue.centres(index,2);
                walker.z = tissue.centres(index,3);
                walker.lab = 1;
                
            elseif(chance >= tissue.ficvf)                
                walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
                walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
                walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
                walker.lab = NaN;
                walker = isintra(walker,tissue);
                while(walker.lab~=0)
                    walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
                    walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
                    walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
                    walker = isintra(walker,tissue);
                end
            end
                
            
        otherwise
            warning('walker initialisation not recognised, defaulting to "random_quartile"');
            walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
            walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
            walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
    end
            
        walker.s = 0.0;        
        walker = scat(walker,tissue); %Calculate directional cosines
        walker.t = max_t;
        walker.alive = 1;
    end

%Compute new direction

    function [walker] = scat(walker,tissue)
        walker.ux = -1 + (2*rand);
        walker.uy = -1 + (2*rand);
        walker.uz = -1 + (2*rand);
        temp = abs(walker.ux) + abs(walker.uy) + abs(walker.uz);
        walker.ux = walker.ux/temp;
        walker.uy = walker.uy/temp;
        walker.uz = walker.uz/temp;
        walker.ux = sqrt(abs(walker.ux)) * sign(walker.ux);
        walker.uy = sqrt(abs(walker.uy)) * sign(walker.uy);
        walker.uz = sqrt(abs(walker.uz)) * sign(walker.uz);
        
        %Need separate scatter function to handle walkers in vessels - It's FLOW time...
        
        if(walker.lab == 2)  %Will only catch after walker label assigned?
            %Calculate the branch that the walker is in
            
            %Nearest node
            [nearest_node, ~] = dsearchn(tissue.nodes, tissue.node_T, [walker.x,walker.y,walker.z]);
            
            %Nodes that bifurcate from the nearest node
            dest_nodes = tissue.vtree.dA(:,nearest_node);
            dest_nodes_coords = tissue.nodes(dest_nodes==1,:);
            d_dists = pdist2(dest_nodes_coords,[walker.x,walker.y,walker.z], 'euclidean');
            d_index = find(dest_nodes);
            dest_nodes = [dest_nodes_coords, d_dists, d_index];
            
            
            
            %Nodes that bifurcate to nearest node
            source_nodes = tissue.vtree.dA(nearest_node,:);
            source_nodes_coords = tissue.nodes(source_nodes==1,:);
            s_dists = pdist2(source_nodes_coords, [walker.x,walker.y,walker.z], 'euclidean');
            s_index = find(source_nodes);
            source_nodes = [source_nodes_coords, s_dists s_index];
            
            
            %Shortest distance of either gives branch?
            if(isempty(dest_nodes))
                nodes = source_nodes;
            elseif(isempty(source_nodes))
                nodes = dest_nodes;
            else
                nodes = [dest_nodes;source_nodes];
            end
            [~,order] = sort(nodes(:,4),'ascend');
            nodes = nodes(order,:);
            nearest_node2 = nodes(1,5); %Second-nearest node
            
            %Use second-nearest node to assign random branch if path is bifurcating.
            
            if(length(find(dest_nodes(:,5)==nearest_node2)) == 1)
                %Pick a destination at random
                nearest_node2 = datasample(dest_nodes(:,5),1);
            end
            
            %Calculate direction of flow using branching order of two nodes
            
            if(tissue.vroute(nearest_node) < tissue.vroute(nearest_node2))
                %Collect cartesian coords of two nodes
                x1 = tissue.vtree.X(nearest_node);
                x2 = tissue.vtree.X(nearest_node2);
                y1 = tissue.vtree.Y(nearest_node);
                y2 = tissue.vtree.Y(nearest_node2);
                z1 = tissue.vtree.Z(nearest_node);
                z2 = tissue.vtree.Z(nearest_node2);
                
                %Calculate directional cosines
                delta_x = x2-x1;
                delta_y = y2-y1;    %Changes in x,y,z
                delta_z = z2-z1;
                
                mag_u = sqrt(delta_x^2 + delta_y^2 + delta_z^2);    %Vector magnitude
                
                walker.ux = delta_x/mag_u;
                walker.uy = delta_y/mag_u;
                walker.uz = delta_z/mag_u;               
           
            else
                %Collect cartesian coords of two nodes
                x1 = tissue.vtree.X(nearest_node2);
                x2 = tissue.vtree.X(nearest_node);
                y1 = tissue.vtree.Y(nearest_node2);
                y2 = tissue.vtree.Y(nearest_node);
                z1 = tissue.vtree.Z(nearest_node2);
                z2 = tissue.vtree.Z(nearest_node);
                
                %Calculate directional cosines
                delta_x = x2-x1;
                delta_y = y2-y1;    %Changes in x,y,z
                delta_z = z2-z1;
                
                mag_u = sqrt(delta_x^2 + delta_y^2 + delta_z^2);    %Vector magnitude
                
                walker.ux = delta_x/mag_u;
                walker.uy = delta_y/mag_u;
                walker.uz = delta_z/mag_u;
            end
            %Add some noise to prevent sticking
            walker.ux = walker.ux + random('norm',0,0.05);
            walker.uy = walker.uy + random('norm',0,0.05); %BJ - Increase vasc diffusivity?
            walker.uz = walker.uz + random('norm',0,0.05);
            %Check to see if walker is at terminal point
            if(isempty(dest_nodes_coords))
                walker.ux = 0;
                walker.uy = 0;
                walker.uz = 0;
            end                            
        end
    end

%Move the walker by stepsize s

    function [walker] = move(walker)
        
        newpos.x = walker.x + walker.sx;
        newpos.y = walker.y + walker.sy;
        newpos.z = walker.z + walker.sz;     
        
        walker.x = newpos.x;
        walker.y = newpos.y;
        walker.z = newpos.z;
        
        walker.t = walker.t - dt;
        
        if(walker.t <= 0)
            walker.alive = 0;
        end
    end

%Implement periodic boundary if walker hits boundary

    function [walker] = voxelboundary(walker)
        
        %X-direction
        if(walker.x >= tissue.dimx)
            walker.x = 0;
        elseif (walker.x <= 0)
            walker.x = tissue.dimx;
        end
        
        %Y-direction
        if(walker.y >= tissue.dimy)
            walker.y = 0;
        elseif(walker.y <= 0)
            walker.y = tissue.dimy;
        end
        
        %Z-direction
        if(walker.z >= tissue.dimz)
            walker.z = 0;
        elseif(walker.z <= 0)
            walker.z = tissue.dimz;
        end
    end

%Adjust step size for boundaries

    function [walker,bhit,thit,vhit] = willhit(walker,tissue)
        bhit = 0;
        thit = 0;
        vhit = 0;
        %Calculate directional step sizes
        walker.s = sqrt(6*tissue.di*dt); %Always 3D
        %----------------------------------------------------------------------
        % If walker is in vessels, need to calculate flow velocity
        %----------------------------------------------------------------------
            if(walker.lab == 2) 
                walker.s = 3e-4 * dt; %3e-4 is typical flow speed in capillaries, according to wikipedia...
            end
        %----------------------------------------------------------------------
        walker.sx = walker.s*walker.ux;
        walker.sy = walker.s*walker.uy;
        walker.sz = walker.s*walker.uz;
        
        %Calculate distance to voxel boundary
        if(walker.ux > 0)
            dtb.x = (tissue.dimx - walker.x)/walker.ux;
        else
            dtb.x = walker.x/walker.ux;
        end
        
        if(walker.uy > 0)
            dtb.y = (tissue.dimy - walker.y)/walker.uy;
        else
            dtb.y = walker.y/walker.uy;
        end
        
        if(walker.uz > 0)
            dtb.z = (tissue.dimz - walker.z)/walker.uz;
        else
            dtb.z = walker.z/walker.uz;
        end
        
        %Fiddle with dtb to ensure it does not =0.
        
        if(dtb.x == 0 || dtb.y == 0 || dtb.z == 0)
            if(dtb.x == 0)
                if(walker.x == 0 && walker.ux < 0)
                    dtb.x = dtb.x-1e-10;
                elseif(walker.x == tissue.dimx && walker.ux > 0)
                    dtb.x = dtb.x+1e-10;
                end
            end
            if(dtb.y == 0)
                if(walker.y == 0 && walker.uy < 0)
                    dtb.y = dtb.y-1e-10;
                elseif(walker.y == tissue.dimy && walker.uy > 0)
                    dtb.y = dtb.y+1e-10;
                end
            end
            if(dtb.z == 0)
                if(walker.z == 0 && walker.uz < 0)
                    dtb.z = dtb.z-1e-10;
                elseif(walker.z == tissue.dimz && walker.uz > 0)
                    dtb.z = dtb.z+1e-10;
                end
            end
        end
        
        next_x = walker.x + walker.sx;
        next_y = walker.y + walker.sy;
        next_z = walker.z + walker.sz;
        
        if(next_x<0 ||next_x>tissue.dimx)
            walker.sx = dtb.x;
            bhit = 1;
        end
        if(next_y<0 || next_y>tissue.dimy)
            walker.sy = dtb.y;
            bhit = 1;
        end
        if(next_z < 0 || next_z>tissue.dimz)
            walker.sz = dtb.z;
            bhit = 1;
        end
        
        if(walker.lab == 1)
            centre = tissue.centres(walker.cellid,:);
            next_pos = [next_x,next_y,next_z];
            dist = pdist([centre;next_pos], 'euclidean');
            if (dist > tissue.rads(walker.cellid))
                thit = 1;
            end
        elseif(walker.lab == 2)
%             if(inpolyhedron(tissue.vhull, [next_x, next_y, next_z]) == 0)
%                 vhit = 1;
%             end
            
        end
    end

%Complete one full step

    function [walker] = step(walker)
        
        walker.didmove = 0;
                    
        %%% CHANGE -  Use step location to decide step behaviour%%%
        
        [walker,bhit,thit,vhit] = willhit(walker,tissue);
        
        if(bhit==1||thit==1||vhit==1)
            if(bhit==1&&vhit==0&&thit==0)
                walker = move(walker);
                walker = voxelboundary(walker);
                walker = scat(walker,tissue); %BJ edit 13/06/16
                walker.didmove = 1; % ????? could this affect diffusion signal?
            elseif(thit==1||vhit==1)
                walker = scat(walker,tissue);
                walker.didmove = 0;
            end
        else
            walker = move(walker);
            walker = scat(walker, tissue);
            walker.didmove = 1;
        end
        
    end

    %Find whether walker is in intracellular, extracellular, or vascular space
    function [walker] = isintra(walker,tissue)
        [k,d] = dsearchn(tissue.centres, tissue.T, [walker.x,walker.y,walker.z]);
        if(walker.lab == 2)
            walker.lab = 2; %vascular
            walker.cellid = NaN;
        elseif(d <= tissue.rads(k))
            walker.lab = 1; %intracellular
            walker.cellid = k;
        else
            walker.lab = 0; %extracellular
            walker.cellid = NaN;
        end
    end

%-----------------------------------------------------------------------------
% Propagate one walker
%-----------------------------------------------------------------------------
tic
h = waitbar(0,['Propagated 0 out of', num2str(numwalkers), ' walkers']);
stats.ficvf = 0;
stats.fee = 0;  %Counters to assess true volume fractions
stats.fvasc = 0;
for i = 1:numwalkers
    walker = setwalker(stats, initialisation);
    walker = isintra(walker,tissue);
    if(walker.lab == 0)
        stats.fee = stats.fee + 1;
    elseif(walker.lab == 1)
        stats.ficvf = stats.ficvf + 1;
    else
        stats.fvasc = stats.fvasc + 1;
    end
    j=1;
    X(i,j) = walker.x;
    Y(i,j) = walker.y;
    Z(i,j) = walker.z;
    S(i,j) = walker.s;
    T(i,j) = max_t - walker.t;
    I(i,j) = i-1;
    j = j+1;
    while (walker.alive == 1)
        walker = step(walker);
        if(walker.didmove == 1)
            X(i,j) = walker.x;
            Y(i,j) = walker.y;
            Z(i,j) = walker.z;
            S(i,j) = walker.s;
            T(i,j) = max_t - walker.t;
            I(i,j) = i-1;
            j = j+1;     
                    
        end
    end
    waitbar(i/numwalkers, h, ['Propagated ', num2str(i),' out of ', num2str(numwalkers), ' walkers']);
    traj(:,:,i) = [T(i,:)',I(i,:)',X(i,:)',Y(i,:)',Z(i,:)',S(i,:)'];  
end
close(h);
%Create headerinfo
headerinfo.duration = max_t;
headerinfo.nwalkers = numwalkers;
headerinfo.tmax = (max_t/dt);
runtime = toc
stats.runtime = runtime;
end