function [traj, headerinfo, stats] = tissuediffusion_amira(tissue, numwalkers, max_t, dt, initialisation)
% Monte Carlo simulation of diffusion MR acquisition from water molecule diffusion through tissue
% using a real vascular network structure.
%
% Input: tissue - A structure containing tissue properties:
%  
%             tissue.dimx     -   Voxel X-dimension
%             tissue.dimy     -   Voxel Y-dimension
%             tissue.dimz     -   Voxel Z-dimension
%             tissue.ficvf    -   Intracellular Volume Fraction
%             tissue.fee      -   Extracellular Volume Fraction
%             tissue.fvasc    -   Vascular Volume Fraction
%             tissue.rads     -   Column-vector of cell radii
%             tissue.voxvol   -   Total Volume of Voxel
%             tissue.cellvol  -   Total Volume of all cells
%             tissue.numcells -   Number of cells in voxel
%             tissue.di       -   Diffusivity of intracellular/extracellular space
%             tissue.dv       -   Diffusivity of vascular space
%             tissue.centres  -   N*3 matrix of cell positions
%             tissue.T        -   Delaunay triangulation of cell locations
%             tissue.cellrad  -   Vector with mean radius and std deviation of cell size
%             tissue.adata    -   Structure containing the segmented vascular network.
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
%             -> 'free_quartile'      - Only initialise in EES, inner half.
% 
% Outputs: traj - Matrix of walker trajectories
%          headerinfo - Structure containing information needed to create trajectory files
%
% Author: Ben Jordan - rmapbjo@ucl.ac.uk

%---------------------------------------------------------------------------------------------------
% Funtions for walker propagation
%---------------------------------------------------------------------------------------------------

%Initialise walker properties
    function [walker, stats] = setwalker(stats, initialisation)
        
        %Initialise walker positions based on user input 'initialisation'
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
                    walker = isintra(walker, tissue);
                    while(walker.lab~=0)
                        walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
                        walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
                        walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
                        walker = isintra(walker, tissue);
                    end
                else
                    index = ceil(rand * length(tissue.adata(1).Val));
                    walker.x = tissue.adata(1).Val(index, 1);
                    walker.y = tissue.adata(1).Val(index, 2);
                    walker.z = tissue.adata(1).Val(index, 3);
                    walker.cvx = index;
                    walker.lab = 2;
                end
                
            case 'random_cells'
                %Initialise only in intercellular space, random cell choice
                index = ceil(rand * length(tissue.centres));
                walker.x = tissue.centres(index, 1);
                walker.y = tissue.centres(index, 2);
                walker.z = tissue.centres(index, 3);
                
            case 'fixed_cells'
                %Initialise only in intercellular space, systematic cell choice
                walker.x = tissue.centres(i,1);
                walker.y = tissue.centres(i,2);
                walker.z = tissue.centres(i,3);
                
            case 'vasc'
                %Initialise only in the vasculature
                index = ceil(rand * length(tissue.adata(1).Val));
                walker.x = tissue.adata(1).Val(index, 1);
                walker.y = tissue.adata(1).Val(index, 2);
                walker.z = tissue.adata(1).Val(index, 3);
                walker.cvx = index;
                walker.lab = 2;
                
            otherwise
                warning('Walker initialisation not recognised, defaulting to "random_quartile"');
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
                    walker = isintra(walker, tissue);
                    while(walker.lab~=0)
                        walker.x = (rand * (tissue.dimx/2)) + (tissue.dimx/4);
                        walker.y = (rand * (tissue.dimy/2)) + (tissue.dimy/4);
                        walker.z = (rand * (tissue.dimz/2)) + (tissue.dimz/4);
                        walker = isintra(walker, tissue);
                    end
                else
                    index = ceil(rand * length(tissue.adata(1).Val));
                    walker.x = tissue.adata(1).Val(index, 1);
                    walker.y = tissue.adata(1).Val(index, 2);
                    walker.z = tissue.adata(1).Val(index, 3);
                    walker.cvx = index; %Current walker vertex
                    walker.lab = 2;
                end
        end
        
        walker.s = 0.0;
        walker = scat(walker, tissue); %Calculate directional cosines
        walker.t = max_t';
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
                       
            %Locate nodes adjacent to current node
            
            nextVertexIdx = find(ismember(tissue.adata(2).Val(:,1), walker.cvx-1, 'rows')); %Nodes are listed from 0, in initial list, and 1 everywhere else for some STUPID reason.
            nextVertex = tissue.adata(2).Val(find(ismember(tissue.adata(2).Val(:,1), walker.cvx-1, 'rows')),2);
            prevVertexIdx = find(ismember(tissue.adata(2).Val(:,2), walker.cvx-1, 'rows'));
            prevVertex = tissue.adata(2).Val(find(ismember(tissue.adata(2).Val(:,2), walker.cvx-1, 'rows')),1);
            %Index of current vertex and adjacent vertices in edge list 
            edgeIndex = sum(tissue.adata(3).Val(1:(find(ismember(tissue.adata(2).Val(:,1), walker.cvx-1, 'rows'),1))))-1; %Not even sure why this -1 is needed... But IT IS.
            
            %Pre-allocate for speed
            nxtVertexEdgeIdx = zeros(size(nextVertex));
            prvVertexEdgeIdx = zeros(size(prevVertex));
            
            for v = 1:length(nextVertex)
                nxtVertexEdgeIdx(v) = sum(tissue.adata(3).Val(1:nextVertexIdx(v)));
            end
            
            for v = 1:length(prevVertex)
                prvVertexEdgeIdx(v) = sum(tissue.adata(3).Val(1:prevVertexIdx(v)))-1; %NOW I KNOW... because it's the source node, and so is indexed first;
            end
            
            %Change - to handle previous vertex also, may flow backwards if
            %pressure dictates.

            if(isempty(nextVertex)) %Node is termination point
                walker.lab = NaN; %Leak out into tissue.
                return;
            elseif(length(nextVertex > 1)) %Node is a bifurcation/trifurcation point

                %Decide which node to propagate to based on pressure drop
                currentVertexPressure = tissue.adata(7).Val(edgeIndex);
                
                nxtVertexPressures = zeros(size(nextVertex));
                prvVertexPressures = zeros(size(prevVertex));
                
                %Pressure drops
                nxtVertexDrops = zeros(size(nextVertex), 2);
                prvVertexDrops = zeros(size(prevVertex), 2);
                
                for v = 1:length(nxtVertexPressures)
                    nxtVertexPressures(v) = tissue.adata(7).Val(nxtVertexEdgeIdx(v));
                    nxtVertexDrops(v, 1) = currentVertexPressure - nxtVertexPressures(v);
                    nxtVertexDrops(v, 2) = nxtVertexEdgeIdx(v);
                end
                
                for v = 1:length(prvVertexPressures)
                    prvVertexPressures(v) = tissue.adata(7).Val(prvVertexEdgeIdx(v));
                    prvVertexDrops(v, 1) = currentVertexPressure - prvVertexPressures(v);
                    prvVertexDrops(v, 2) = prvVertexEdgeIdx(v);
                end
                               
                
                pDrops = [nxtVertexDrops; prvVertexDrops];%Concatenate pressure drops
                [~,ind] = min(pDrops(:,1));
                destIdx = pDrops(ind,2);
                
            else %Node is continuous
                prevVertexPressure = tissue.adata(7).Val(prvVertexEdgeIdx);
                currentVertexPressure = tissue.adata(7).Val(edgeIndex);
                nextVertexPressure = tissue.adata(7).Val(edgeIndex + tissue.adata(3).Val(nextVertexIdx));
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
%Create function handles for parallel loop
set = @setwalker;
ii = @isintra;
stp = @step;

stats.runtime = 0;

%Initialise outputs

X = zeros(numwalkers, max_t/dt+2);
Y = zeros(numwalkers, max_t/dt+2);
Z = zeros(numwalkers, max_t/dt+2);
S = zeros(numwalkers, max_t/dt+2);
T = zeros(numwalkers, max_t/dt+2);
I = zeros(numwalkers, max_t/dt+2);

parfor i = 1:numwalkers
    
    Xtemp = zeros(1,max_t/dt);
    Ytemp = zeros(1,max_t/dt);
    Ztemp = zeros(1,max_t/dt);
    Stemp = zeros(1,max_t/dt);
    Ttemp = zeros(1,max_t/dt);
    Itemp = zeros(1,max_t/dt);
    walker = set(stats, initialisation);
    walker = ii(walker,tissue);
    j=1;
    Xtemp(j) = walker.x;
    Ytemp(j) = walker.y;
    Ztemp(j) = walker.z;
    Stemp(j) = walker.s;
    Ttemp(j) = max_t - walker.t;
    Itemp(j) = i-1;
    j = j+1;
    while (walker.alive == 1)
        walker = stp(walker);
        if(walker.didmove == 1)
            Xtemp(j) = walker.x;
            Ytemp(j) = walker.y;
            Ztemp(j) = walker.z;
            Stemp(j) = walker.s;
            Ttemp(j) = max_t - walker.t;
            Itemp(j) = i-1;
            j = j+1;
            
        end
    end
    X(i,:) = Xtemp(:);
    Y(i,:) = Ytemp(:);
    Z(i,:) = Ztemp(:);
    S(i,:) = Stemp(:);
    T(i,:) = Ttemp(:);
    I(i,:) = Itemp(:);
    
end

for i = 1:numwalkers
    traj(:,:,i) = [T(i,:)',I(i,:)',X(i,:)',Y(i,:)',Z(i,:)',S(i,:)'];
end
%close(h);
%Create headerinfo
headerinfo.duration = max_t;
headerinfo.nwalkers = numwalkers;
headerinfo.tmax = (max_t/dt);
runtime = toc
stats.runtime = runtime;
end