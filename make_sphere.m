% make a sphere mesh


[V,Tri,~,Ue]=ParticleSampleSphere('N',30); 
R = 12E-6;
V_size = V.*R;
for i = 1:length(V)
    Data.vertex.x(i) = V_size(i,1);
    Data.vertex.y(i) = V_size(i,2);
    Data.vertex.z(i) = V_size(i,3);
end

for i = 1:length(Tri)
    Data.face.vertex_indices{i} = Tri(i,:)-1;  
end
figure, h=trimesh(Tri,V(:,1),V(:,2),V(:,3)); set(h,'EdgeColor','b'), axis equal
ply_write(Data,['sphere_rad' num2str(R*1E6) '.ply'],'ascii');
