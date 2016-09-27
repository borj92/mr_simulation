function [p, t] = TriSphere(N, R)
% TRISPHERE: Returns the triangulated model of a sphere using the
% icosaedron subdivision method.
%
% INPUT:
% N (integer number) indicates the number of subdivisions,
%   it can assumes values between 0-inf. The greater N the better will look
%   the surface but the more time will be spent in surface costruction and
%   more triangles will be put in the output model.
%
% OUTPUT:
% In p (nx3) and t(mx3) we can find points and triangles indexes
% of the model. The sphere is supposed to be of unit radius and centered in
% (0,0,0). To obtain spheres centered in different location, or with
% different radius, is just necessary a traslation and a scaling
% trasformation. These operation are not perfomed by this code beacuse it is
% extrimely convinient, in order of time perfomances, to do this operation
% out of the function avoiding to call the costruction step each time.
%
% NOTE:
% This function is more efficient than the matlab command sphere in
% terms of dimension fo the model/ accuracy of recostruction. This due to
% well traingulated model that requires a minor number of patches for the
% same geometrical recostruction accuracy. Possible improvement are possible
% in time execution and model subdivision flexibilty.
%
% EXAMPLE:
%
%  N=5;
%
%  [p,t] = TriSphere(N);
%
%  figure(1) axis equal hold on trisurf(t,p(:,1),p(:,2),p(:,3)); axis vis3d
%  view(3)

% Author: Giaccari Luigi Created:25/04/2009%
% For info/bugs/questions/suggestions : giaccariluigi@msn.com
% ORIGINAL NAME: BUILDSPHERE
%
% Adjusted by Rody Oldenhuis (speed/readability)

    % error traps
    error(nargoutchk(1,2,nargout));
    if ~isscalar(N)
        error('Buildsphere:N_mustbe_scalar',...
            'Input N must be a scalar.');
    end    
    if round(N) ~= N
        error('Buildsphere:N_mustbe_scalar',...
            'Input N must be an integer value.');
    end

    % Coordinates have been taken from Jon Leech' files

    % Twelve vertices of icosahedron on unit sphere
    tau = 0.8506508083520400; % t   = (1+sqrt(5))/2, tau = t/sqrt(1+t^2)
    one = 0.5257311121191336; % one = 1/sqrt(1+t^2)  (unit sphere)    
    p = [
        +tau  +one  +0     % ZA
        -tau  +one  +0     % ZB
        -tau  -one  +0     % ZC
        +tau  -one  +0     % ZD
        +one  +0    +tau   % YA
        +one  +0    -tau   % YB
        -one  +0    -tau   % YC
        -one  +0    +tau   % YD
        +0    +tau  +one   % XA
        +0    -tau  +one   % XB
        +0    -tau  -one   % XC
        +0    +tau  -one]; % XD

    % Structure for unit icosahedron
    t = [  
         5  8  9 
         5 10  8 
         6 12  7 
         6  7 11 
         1  4  5 
         1  6  4 
         3  2  8 
         3  7  2 
         9 12  1 
         9  2 12 
        10  4 11 
        10 11  3 
         9  1  5 
        12  6  1 
         5  4 10 
         6 11  4 
         8  2  9 
         7 12  2 
         8 10  3 
         7  3 11 ];

    % possible quick exit
    if N == 0, return, end

    % load pre-generated trispheres (up to 8 now...)
    if N <= 8
        S = load(['TriSphere', num2str(N), '.mat'],'pts','idx');
        p = S.pts; t = S.idx; 
        if nargin == 2, p = p*R; end
        return
    else
        % if even more is requested (why on Earth would you?!), make sure to START 
        % from the maximum pre-loadable trisphere
        S = load('TriSphere8.mat','pts','idx');
        p = S.pts; t = S.idx; clear S; N0 = 10;
    end

    % how many triangles/vertices do we have? 
    nt = size(t,1); np = size(p,1); totp = np;    
    % calculate the final number of points    
    for ii=N0:N, totp = 4*totp - 6; end    
    % initialize points array
    p = [p; zeros(totp-12, 3)];

    % determine the appropriate class for the triangulation indices
    numbits   = 2^ceil(log(log(totp+1)/log(2))/log(2));
    castToInt = ['uint',num2str(numbits)];

    % issue warning when required
    if numbits > 64
        warning('TriSphere:too_many_notes',...
            ['Given number of iterations would require a %s class to accurately ',...
            'represent the triangulation indices. Using double instead; Expect ',...
            'strange results!']);
        castToInt = @double;
    else
        castToInt = str2func(castToInt);
    end

    % refine icosahedron N times
    for ii = N0:N
        % initialize inner loop
        told  = t;
        t = zeros(nt*4, 3);
        % Use sparse. Yes, its slower in a loop, but for N = 6 the size is
        % already ~10,000x10,000, growing by a factor of 4 with every
        % increasing N; its simply too memory intensive to use zeros().
        peMap = sparse(np,np); 
        ct    = 1;        
        % loop trough all old triangles        
        for j = 1:nt

            % some helper variables
            p1 = told(j,1);
            p2 = told(j,2);
            p3 = told(j,3);
            x1 = p(p1,1); x2 = p(p2,1); x3 = p(p3,1);
            y1 = p(p1,2); y2 = p(p2,2); y3 = p(p3,2);
            z1 = p(p1,3); z2 = p(p2,3); z3 = p(p3,3);

            % first edge
            % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            % preserve triangle orientation
            if p1 < p2, p1m = p1; p2m = p2; else p2m = p1; p1m = p2; end

            % If the point does not exist yet, calculate the new point            
            p4 = peMap(p1m,p2m);
            if p4 == 0
                np = np+1;
                p4 = np;
                peMap(p1m,p2m) = np;%#ok
                p(np,1) = (x1+x2)/2;
                p(np,2) = (y1+y2)/2; 
                p(np,3) = (z1+z2)/2;
            end

            % second edge
            % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            % preserve triangle orientation
            if p2 < p3; p2m = p2; p3m = p3; else p2m = p3; p3m = p2; end

            % If the point does not exist yet, calculate the new point
            p5 = peMap(p2m,p3m);
            if p5 == 0
                np = np+1;
                p5 = np;
                peMap(p2m,p3m) = np;%#ok
                p(np,1) = (x2+x3)/2; 
                p(np,2) = (y2+y3)/2; 
                p(np,3) = (z2+z3)/2;
            end

            % third edge
            % -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

            % preserve triangle orientation
            if p1 < p3; p1m = p1; p3m = p3; else p3m = p1; p1m = p3; end

            % If the point does not exist yet, calculate the new point            
            p6 = peMap(p1m,p3m);
            if p6 == 0
                np = np+1;
                p6 = np;
                peMap(p1m,p3m) = np;%#ok
                p(np,1) = (x1+x3)/2; 
                p(np,2) = (y1+y3)/2;
                p(np,3) = (z1+z3)/2;                
            end

            % allocate new triangles
            %   refine indexing
            %          p1
            %          /\
            %         /t1\
            %      p6/____\p4
            %       /\    /\
            %      /t4\t2/t3\
            %     /____\/____\
            %    p3    p5     p2            
            t(ct,1) = p1; t(ct,2) = p4; t(ct,3) = p6; ct = ct+1;            
            t(ct,1) = p4; t(ct,2) = p5; t(ct,3) = p6; ct = ct+1;
            t(ct,1) = p4; t(ct,2) = p2; t(ct,3) = p5; ct = ct+1;            
            t(ct,1) = p6; t(ct,2) = p5; t(ct,3) = p3; ct = ct+1;

        end % end subloop
        % update number of triangles
        nt = ct-1;         
    end % end main loop

    % normalize all points to 1 (or R)
    p = bsxfun(@rdivide, p, sqrt(sum(p.^2,2)));
    if (nargin == 2), p = p*R; end
    % convert t to proper integer class
    t = castToInt(t);    

end % funciton TriSphere