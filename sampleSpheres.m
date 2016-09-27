function [centers, rads] = sampleSpheres( dims, n, radius) 
 % main function which is to be called for adding multiple spheres in a cube
 % dims is assumed to be a row vector of size 1-by-ndim


 % preallocate
 ndim = numel(dims);
 centers = zeros( n, ndim );
 rads = zeros( n, 1 );
 ii = 1;
 %h = waitbar(0, 'Populating Spheres');
 while ii <= n
      [centers(ii,:), rads(ii) ] = randomSphere(dims, radius);     
      if nonOver( centers(1:ii,:), rads(1:ii) )
           ii = ii + 1; % accept and move on
      end
      %waitbar(ii/n);
%       prog = [num2str(ii), '/', num2str(n)];
%       prog
 end
 %close(h)
 
end
