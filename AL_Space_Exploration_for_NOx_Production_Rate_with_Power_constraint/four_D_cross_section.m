load('matlab_BO.mat')
his_list_y = his_list_y.*(60*24*1e-3/(14*2.5*1e-6));


O2_lin = linspace(0,0.2095,250);
dis_lin = linspace(2,15,250);
[O2frac, distance] = meshgrid(O2_lin, dis_lin);

flat_point = [repmat(his_list_x(93,1),size(reshape(O2frac,[],1))), reshape(O2frac,[],1), reshape(distance,[],1), repmat(his_list_x(93,4),size(reshape(O2frac,[],1)))];

enr_flat = ftir_NN_4D(flat_point).*(60*24*1e-3/(14*2.5*1e-6));

out_enr = reshape(enr_flat,[],length(O2_lin));

[max_val,max_id] = max(1./enr_flat);

hold on
box on
surf(O2frac,distance,1./out_enr)
mybar = colorbar;
colormap(hot) 

scatter3(O2frac(max_id),distance(max_id),max_val,500,'k','filled','p')
%scatter3(his_list_x(1:current_step,2),his_list_x(1:current_step,3),his_list_y(1:current_step),50,'r','filled','o')
%scatter3(his_list_x(current_step+1,2),his_list_x(current_step+1,3),his_list_y(current_step+1),500,'k','filled','p')

view(2)

shading interp
xlim([0,0.2095])
ylim([2,15])
xlabel('O_2 mole fraction (%)','fontsize',14)
ylabel('Distance(mm)','fontsize',14)
mybar.Label.String = '1/E_{N_r}(tN/GJ)';
set(gca,'Fontsize',14)
mybar.Label.FontSize = 18;

set(gcf,'renderer','Painters')
