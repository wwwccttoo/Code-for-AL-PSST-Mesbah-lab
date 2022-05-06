load('for_three_d_plot_again.mat');
dt = delaunayTriangulation(mm',mA') ;
tri = dt.ConnectivityList;
xi = dt.Points(:,1);
yi = dt.Points(:,2);

F = scatteredInterpolant(mm',mA',energy_cost');
zi = F(xi,yi) ;
trisurf(tri,xi,yi,zi,'FaceColor','r','FaceAlpha',0.5) 
xlabel('Distance (mm)');
ylabel('Current (mA)');
zlabel('Energy cost (GJ/tN)');
colorbar
view(2)
shading interp