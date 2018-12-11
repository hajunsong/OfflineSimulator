clc; clear all;

map = load('../map_10cm_anim.csv');
path = load('../3d_path_10cm.csv');

wp = size(east,2);
row = size(map,1);
col = size(map,2);

gap = 0.1;
x = linspace(0,1080,col);
y = linspace(0,1170,row);
[X,Y] = meshgrid(x,y);

figure
set(gcf,'Color',[1,1,1])
% imagesc(x,y,map);
mesh(X,Y,map);
hold on
% plot(east*gap,north*gap,'r','LineWidth',2);
plot3(path(:,1), path(:,2), path(:,3),'r','LineWidth',2)
index = 1270;
plot3(path(index,1), path(index,2), path(index,3),'ko','LineWidth',2,'MarkerSize',5,'MarkerFaceColor','k')
grid on
range = 5;
axis([path(index,1)-range,path(index,1)+range, path(index,2)-range,path(index,2)+range ,path(index,3)-7,path(index,3)+7])
xlabel('X')
ylabel('Y')
set(gca,'FontSize',13)
% view(160,50)