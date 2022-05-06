x = linspace(-2,10,100);
y = demof_for_GP(x);%+normrnd(0,0.05,[1,100]);

load('GP_train_demo.mat');
load('previous11evals.mat');

current_xtrain = previous_11_iters;
current_ytrain = demof_for_GP(current_xtrain);


uqlab
% train GP
MetaOpts.Type = 'Metamodel';
MetaOpts.Scaling=1;
MetaOpts.MetaType = 'Kriging';
MetaOpts.Corr.Family = 'Gaussian';
MetaOpts.Corr.Isotropic = 1;

%MetaOpts.ExpDesign.X = [current_xtrain([2,5,7,10,11]);9;9.8;1.8;0.9;2.402040204020402;3.828982898289829];
%MetaOpts.ExpDesign.Y = [current_ytrain([2,5,7,10,11]);demof_for_GP([9;9.8;1.8;0.9;2.402040204020402;3.828982898289829])];

MetaOpts.ExpDesign.X = [current_xtrain([2,5,7,4]);0.9;9.9;9.0483];
MetaOpts.ExpDesign.Y = [current_ytrain([2,5,7,4]);demof_for_GP([0.9;9.9;9.0483])];
myKrigingMat = uq_createModel(MetaOpts);


expected_improvement_cal(10,myKrigingMat,current_ytrain)

EIfunc = @(x)expected_improvement_cal(x,myKrigingMat,current_ytrain)*-1;



[Ypred,Yvar] = uq_evalModel(myKrigingMat,x');

Lwidth = 2;

t = tiledlayout(3,3);
box on

ax1 = nexttile(1,[2 3]);
xconf = [x x(end:-1:1)] ;         
yconf = [(Ypred+2*sqrt(Yvar))' (Ypred(end:-1:1)-2*sqrt(Yvar(end:-1:1)))'];
h1 = fill(ax1,xconf,yconf,'red','EdgeColor','None');
h1.FaceColor = [0,0.475,0.698];      
h1.FaceAlpha = 0.3;

hold on

ylim([-0.05,1.4])
h2 = plot(ax1,x',y','--','linewidth', Lwidth, 'color','black');

%[trainx_sorted,trainx_sorted_id] = sort(current_xtrain);
h3 = plot(ax1,x',Ypred,'linewidth', Lwidth,'color','blue');
h4 = scatter(ax1,MetaOpts.ExpDesign.X, MetaOpts.ExpDesign.Y,100,'red','.');

%hl = legend('2\sigma confidence interval','Target function','GP prediction','Training data');
%set(hl,'box','off')

ylabel('$f(x)$','fontsize',20,'interpreter','latex','fontweight','bold')
title('8 training points','fontsize',20)
set(gca,'fontweight','bold','Fontsize',12)
hold off


ax2 = nexttile(7,[1 3]);
box on

hold on
h5 = plot(ax2,linspace(-2,10,10000)',-1*EIfunc(linspace(-2,10,10000)'),'linewidth', Lwidth,'color',[0 0.5 0]);
[real_min,min_id]=min(EIfunc(linspace(-2,10,10000)'));
x_inter = linspace(-2,10,10000);
h6 = plot(ax2,[x_inter(min_id) x_inter(min_id)],[0,0.01 ],'--','linewidth', Lwidth,'color',[0.9290, 0.6940, 0.1250]);
h7 = scatter(ax2,x_inter(min_id),0,500,[0.9290, 0.6940, 0.1250],'filled','p');
%set(gca,'fontweight','bold','Fontsize',12)
ylabel('$\alpha_{EI}(x)$','fontsize',20,'interpreter','latex','fontweight','bold')

linkaxes([ax1,ax2],'x');
xticklabels(ax1,{})
xlabel('$x$','fontsize',20,'interpreter','latex','fontweight','bold')

hl = legend([h1,h2,h3,h4,h5,h7],'2\sigma confidence interval','Target function','GP prediction','Training data','Acquisition function value','Suggested next point');
set(hl,'box','off')
set(gca,'fontweight','bold','Fontsize',12)

function EI_value = expected_improvement_cal(x,model,current_y)
    [current_mean,current_variance] = uq_evalModel(model,x);
    current_min = min(current_y);
    EI_value = (current_min-current_mean).*normcdf((current_min-current_mean)./sqrt(current_variance))+sqrt(current_variance).*normpdf((current_min-current_mean)./sqrt(current_variance));
end
