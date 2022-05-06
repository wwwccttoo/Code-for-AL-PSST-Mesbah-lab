Fontsize  = 14;

x_bar = ['I_{dis}', 'O_2%', 'd','Q'];
y_energycost = [0.02700792, 0.57121186, 0.3444228, 0.14375612]; %energy cost
y_NOx = [0.12233159 0.32637085 0.55963974 0.11571165];
y_power = [0.35428996 0.00458087 0.63638527 0.02686357];
b = bar([y_energycost;y_NOx;y_power] ,'facecolor', 'flat');
clr = [0,0.475,0.698;0.173,0.627,0.173;0.749,0.749,0;0,0.749,0.749];
b.CData = clr;

set(gca(),'FontSize',Fontsize,'xticklabel',{'$E_{\mathrm{N}_r}$','$QC_{\mathrm{NO}_\mathrm{x}}$', '$P_\mathrm{dis}$'},'TickLabelInterpreter','latex')
ylabel("$S_T$",'Interpreter','latex' )
lg = legend('$I_{\mathrm{dis}}$', '$\mathrm{O}_2\%$', '$d$','$Q$','Interpreter','latex');
set(lg,'box','off')
