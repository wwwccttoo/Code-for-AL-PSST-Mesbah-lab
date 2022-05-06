x = linspace(-2,10,100);
y = demof_for_GP(x);%+normrnd(0,0.05,[1,100]);

xtrain = [7*rand(1,20)-1.5,rand(1,15)*2.5+7];
ytrain = demof_for_GP(xtrain)+normrnd(0,0.02,[1,35]);


uqlab
% train GP
MetaOpts.Type = 'Metamodel';
MetaOpts.Scaling=1;
MetaOpts.MetaType = 'Kriging';
MetaOpts.Corr.Family = 'Gaussian';
MetaOpts.Corr.Isotropic = 1;

MetaOpts.ExpDesign.X = xtrain';
MetaOpts.ExpDesign.Y = ytrain';
myKrigingMat = uq_createModel(MetaOpts);

[Ypred,Yvar] = uq_evalModel(myKrigingMat,x');

Lwidth = 2;

figure
hold on
box on

xconf = [x x(end:-1:1)] ;         
yconf = [(Ypred+2*sqrt(Yvar))' (Ypred(end:-1:1)-2*sqrt(Yvar(end:-1:1)))'];
h1 = fill(xconf,yconf,'red','EdgeColor','None');
h1.FaceColor = [0,0.475,0.698];      
h1.FaceAlpha = 0.3;
ylim([0,1])
plot(x',y','--','linewidth', Lwidth, 'color','black')
plot(x',Ypred,'linewidth', Lwidth,'color','blue')
scatter(xtrain, ytrain,100,'red','.')


ylabel('$f(x)$','fontsize',20,'interpreter','latex','fontweight','bold')
xlabel('$x$','fontsize',20,'interpreter','latex','fontweight','bold')

hl = legend('2\sigma confidence interval','Target function','GP prediction','Training data');
set(hl,'box','off')
set(gca,'fontweight','bold','Fontsize',12)

