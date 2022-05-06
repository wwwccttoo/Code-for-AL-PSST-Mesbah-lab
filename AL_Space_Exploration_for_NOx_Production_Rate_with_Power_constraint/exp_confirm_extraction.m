clear all
current_distance = openfig('ENr_with_O2frac_distance.fig');
dataObjsY = findobj(current_distance,'-property','YData');

y1 = dataObjsY(1).YData;
y3 = dataObjsY(3).YData;
y4 = dataObjsY(4).YData;
y5 = dataObjsY(5).YData;
y6 = dataObjsY(6).YData;
y2 = dataObjsY(2).YData;

dataObjsX = findobj(current_distance,'-property','XData');
x1 = dataObjsX(1).XData;
x3 = dataObjsX(3).XData;
x4 = dataObjsX(4).XData;
x5 = dataObjsX(5).XData;
x6 = dataObjsX(6).XData;
x2 = dataObjsX(2).XData;
