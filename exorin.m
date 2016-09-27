numruns = 1000;
exorintra = zeros(numruns,1);

for i = 1:numruns
    walker.x = rand * tissue.dimx;
    walker.y = rand * tissue.dimy;
    walker.z = rand * tissue.dimz;
    
    exorintra(i) = inpolyhedron(hull, [walker.x,walker.y,walker.z]);
    %[k,d] = dsearchn(tissue.centres, T, [walker.x,walker.y,walker.z]);
%     if (d <= tissue.rads(k))
%         exorintra(i) = 1;
%     end
end
        