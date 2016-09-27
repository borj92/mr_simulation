function [tissue] = amira2tree(tissue, an)
%% Convert amira-formatted networks into Matlab MC compatible network

tissue.vtree.X = an(1).Val(:,1);
tissue.vtree.Y = an(1).Val(:,2);
tissue.vtree.Z = an(1).Val(:,3);

tissue.vtree.dA = sparse(zeros(length(tissue.vtree.X), length(tissue.vtree.X)));

for i = 1:length(an(2).Val(:,1))
    tissue.vtree.dA(an(2).Val(i,2)+1, an(2).Val(i,1)+1) = 1;
end

tissue.vtree.D = tissue.vtree.D(1:length(tissue.vtree.X));
tissue.vtree.R = tissue.vtree.R(1:length(tissue.vtree.X));

end