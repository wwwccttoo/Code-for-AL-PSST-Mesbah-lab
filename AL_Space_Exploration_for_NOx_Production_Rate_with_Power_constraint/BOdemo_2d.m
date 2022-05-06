load('twoD_demo_starter.mat');

%current_xtrain = [repmat(50,length(points),1),points,repmat(7,length(points),1)];
current_xtrain = points;
current_ytrain = double(evals');


uqlab
% train GP
MetaOpts.Type = 'Metamodel';
MetaOpts.Scaling=1;
MetaOpts.MetaType = 'Kriging';
MetaOpts.Corr.Family = 'Gaussian';
MetaOpts.Corr.Isotropic = 1;

his = 2;
his_list_x = [his_list_x;[O2frac_flat(max_ind),distance_flat(max_ind)]];
his_list_y = [his_list_y;ftir_NN_2D([O2frac_flat(max_ind),distance_flat(max_ind)])];

MetaOpts.ExpDesign.X = his_list_x;
MetaOpts.ExpDesign.Y = his_list_y;
myKrigingMat = uq_createModel(MetaOpts);


%expected_improvement_cal(10,myKrigingMat,current_ytrain)

EIfunc = @(x)expected_improvement_cal(x,myKrigingMat,current_ytrain)*-1;

O2_lin = linspace(0,0.2095,1000);
dis_lin = linspace(2,15,1000);
[O2frac, distance] = meshgrid(O2_lin, dis_lin);

out_enr = zeros(size(O2frac));

out_enr = reshape(-EIfunc([reshape(O2frac,[],1), reshape(distance,[],1)]),[],length(O2_lin));

O2frac_flat = reshape(O2frac,[],1);
distance_flat = reshape(distance,[],1);
out_enr_flat = reshape(out_enr,[],1);
[~,max_ind] = max(out_enr_flat);


hold on
box on
surf(O2frac,distance,out_enr,'FaceAlpha',0.5)
mybar = colorbar;
scatter3(his_list_x(:,1),his_list_x(:,2),his_list_y,50,'r','filled','o')
scatter3(O2frac_flat(max_ind),distance_flat(max_ind),out_enr_flat(max_ind),500,'g','filled','p')

view(2)

shading interp
xlim([0,0.2095])
ylim([2,15])
xlabel('O_2(%)','fontsize',20,'fontweight','bold')
ylabel('Distance(mm)','fontsize',20,'fontweight','bold')
mybar.Label.String = '\alpha_{EI}';
set(gca,'fontweight','bold','Fontsize',12)


