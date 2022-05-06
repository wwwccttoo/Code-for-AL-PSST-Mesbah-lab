load('matlab_BO.mat')

uqlab
% train GP
MetaOpts.Type = 'Metamodel';
MetaOpts.Scaling = 1;
MetaOpts.MetaType = 'Kriging';
MetaOpts.Corr.Family = 'Gaussian';
MetaOpts.Corr.Isotropic = 1;

current_step = 49;

his_list_y = his_list_y.*(60*24*1e-3/(14*2.5*1e-6));

MetaOpts.ExpDesign.X = his_list_x(1:current_step,:);
MetaOpts.ExpDesign.Y = his_list_y(1:current_step,:);
myKrigingMat = uq_createModel(MetaOpts);

EIfunc = @(x)expected_improvement_cal(x, myKrigingMat, MetaOpts.ExpDesign.Y)*-1;

O2_lin = linspace(0,0.2095,300);
dis_lin = linspace(2,15,300);
[O2frac, distance] = meshgrid(O2_lin, dis_lin);

flat_point = [repmat(his_list_x(current_step+1,1),size(reshape(O2frac,[],1))), reshape(O2frac,[],1), reshape(distance,[],1), repmat(his_list_x(current_step+1,4),size(reshape(O2frac,[],1)))];

out_enr = reshape(-EIfunc(flat_point),[],length(O2_lin));


hold on
box on
surf(O2frac,distance,out_enr)
mybar = colorbar;
colormap(jet) 

scatter3(his_list_x(1:current_step,2),his_list_x(1:current_step,3),his_list_y(1:current_step),50,'r','filled','o')
scatter3(his_list_x(current_step+1,2),his_list_x(current_step+1,3),his_list_y(current_step+1),500,'k','filled','p')

view(2)

shading interp
xlim([0,0.2095])
ylim([2,15])
xlabel('O_2 mole fraction (%)','fontsize',14)
ylabel('Distance (mm)','fontsize',14)
mybar.Label.String = '\alpha_{EI}';
set(gca,'Fontsize',14)
mybar.Label.FontSize = 18;

set(gcf,'renderer','Painters')
