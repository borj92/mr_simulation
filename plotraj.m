function output = plotraj(traj,tissue,makevid)
    for i = 1:size(traj,3)
        if(i==1)
            f = figure;
            plot3(traj(:,3,i), traj(:,4,i), traj(:,5,i));
            xlim([0,tissue.dimx]);
            ylim([0,tissue.dimy]);
            zlim([0,tissue.dimz]);
            h = gca;
            box on;
            h.BoxStyle = 'full';
            h.LineWidth = 1;
            h.Color = [0,0,0];
            f.Color = 'none';
            h.XColor = [0.9,0.9,0.9];
            h.YColor = [0.9,0.9,0.9];
            h.ZColor = [0.9,0.9,0.9];
            axis vis3d
            hold on
        else
            plot3(traj(:,3,i), traj(:,4,i), traj(:,5,i))
        end
    end
    hull = hull_tree(tissue.vtree, 8e-6, [],[],[],'-s');
    if(makevid == 1)
        filename = input('Enter filename: ', 's');
        options.FrameRate = 30; options.Duration = 10; options.Periodic = true;
        CaptureFigVid([-20,10;-110,10;-200,10;-290,10;-380,10],filename,options);
    end
end