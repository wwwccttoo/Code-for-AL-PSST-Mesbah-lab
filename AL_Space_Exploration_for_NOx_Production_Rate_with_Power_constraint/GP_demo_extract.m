clear all
BO_RS_store = openfig('GP_demo_plot_light.fig');
dataObjsY = findobj(BO_RS_store,'-property','YData');

y1 = dataObjsY(1).YData;
y2 = dataObjsY(2).YData;
y3 = dataObjsY(3).YData;
y4 = dataObjsY(4).YData;

dataObjsX = findobj(BO_RS_store,'-property','XData');
x1 = dataObjsX(1).XData;
x2 = dataObjsX(2).XData;
x3 = dataObjsX(3).XData;
x4 = dataObjsX(4).XData;
