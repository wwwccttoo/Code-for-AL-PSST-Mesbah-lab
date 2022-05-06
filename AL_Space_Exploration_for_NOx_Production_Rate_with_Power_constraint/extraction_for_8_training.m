BO_demon = openfig('training_7_BO_demo_new.fig');
dataObjsY = findobj(BO_demon,'-property','YData');

y1 = dataObjsY(1).YData;
y2 = dataObjsY(2).YData;
y3 = dataObjsY(3).YData;
y4 = dataObjsY(4).YData;
y5 = dataObjsY(5).YData;
y6 = dataObjsY(6).YData;
y7 = dataObjsY(7).YData;

dataObjsX = findobj(BO_demon,'-property','XData');
x1 = dataObjsX(1).XData;
x2 = dataObjsX(2).XData;
x3 = dataObjsX(3).XData;
x4 = dataObjsX(4).XData;
x5 = dataObjsX(5).XData;
x6 = dataObjsX(6).XData;
x7 = dataObjsX(7).XData;