load('fourD_demo_starter.mat');

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

%his = 2;
%his_list_x = [his_list_x;[O2frac_flat(max_ind),distance_flat(max_ind)]];
%his_list_y = [his_list_y;ftir_NN_2D([O2frac_flat(max_ind),distance_flat(max_ind)])];

his_list_x = current_xtrain(1,:);
his_list_y = current_ytrain(1,:);

rand_second_point = [20+50*rand(1), 0.0005+0.2090*rand(1), 2+13*rand(1), 0.5+6.5*rand(1)];
his_list_x = [his_list_x; rand_second_point];
his_list_y = [his_list_y; ftir_NN_4D(rand_second_point)];



%expected_improvement_cal(10,myKrigingMat,current_ytrain)
for i = 1:150
    tic
    MetaOpts.ExpDesign.X = his_list_x;
    MetaOpts.ExpDesign.Y = his_list_y;
    myKrigingMat = uq_createModel(MetaOpts);
    
    EIfunc = @(x)expected_improvement_cal(x,myKrigingMat,his_list_y)*-1;

    mA_lin = linspace(20,70,15);
    O2_lin = linspace(0,0.2095,15);
    dis_lin = linspace(2,15,15);
    slm_lin = linspace(0.5,7,15);

    [current, O2frac, distance, flowrate] = ndgrid(mA_lin, O2_lin, dis_lin, slm_lin);

    %out_enr = zeros(size(O2frac));

    %out_enr = reshape(-EIfunc([reshape(O2frac,[],1), reshape(distance,[],1)]),[],length(O2_lin));

    %O2frac_flat = reshape(O2frac,[],1);

    out_enr_flat = -EIfunc([reshape(current,[],1), reshape(O2frac,[],1), reshape(distance,[],1), reshape(flowrate,[],1)]);

    currentr_flat = reshape(current,[],1);
    O2frac_flat = reshape(O2frac,[],1);
    distance_flat = reshape(distance,[],1);
    flowrate_flat = reshape(flowrate,[],1);
    [~,max_ind] = max(out_enr_flat);
    
    gs = GlobalSearch;
    problem = createOptimProblem('fmincon','x0',[currentr_flat(max_ind), O2frac_flat(max_ind), distance_flat(max_ind), flowrate_flat(max_ind)],...
        'objective',EIfunc,'lb',[20,0,2,0.5],'ub',[70,0.2095,15,7]);
    find_min = run(gs,problem);
    his_list_x = [his_list_x; find_min];
    his_list_y = [his_list_y; ftir_NN_4D(find_min)];
    toc
end
%{
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

%}
