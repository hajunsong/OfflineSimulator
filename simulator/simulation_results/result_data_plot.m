clc; clear all; close all;

data_indx = 0;
[data, txt] = xlsread(sprintf('%d/SimulationResult.csv', data_indx));

TimeIndex = 1;
VelocityIndex = 29;
WaypointIndex = 30;
DesVelocityIndex = 31;
RoadHeightRFIndex = 37;
RoadHeightRMIndex = 38;
RoadHeightRRIndex = 39;
RoadHeightLFIndex = 40;
RoadHeightLMIndex = 41;
RoadHeightLRIndex = 42;
SuspensionRFIndex = 14;
SuspensionRMIndex = 15;
SuspensionRRIndex = 16;
SuspensionLFIndex = 17;
SuspensionLMIndex = 18;
SuspensionLRIndex = 19;
RF_tire_x = 43;RF_tire_y = 44;RF_tire_z = 45;
RM_tire_x = 46;RM_tire_y = 47;RM_tire_z = 48;
RR_tire_x = 49;RR_tire_y = 50;RR_tire_z = 51;
LF_tire_x = 52;LF_tire_y = 53;LF_tire_z = 54;
LM_tire_x = 55;LM_tire_y = 56;LM_tire_z = 57;
LR_tire_x = 58;LR_tire_y = 59;LR_tire_z = 60;

figure
set(gcf,'Color',[1,1,1])
subplot(211)
plot(data(:,TimeIndex), data(:,VelocityIndex),'LineWidth',2)
grid on
xlabel('Time [s]')
ylabel('Velocity [m/s]')
set(gca,'FontSize',13)

subplot(212)
plot(data(:,WaypointIndex), data(:,VelocityIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Velocity [m/s]')
set(gca,'FontSize',13)

figure
set(gcf,'Color',[1,1,1])
subplot(231)
plot(data(:,WaypointIndex), data(:,RoadHeightRFIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Road Height RF [m]')
set(gca,'FontSize',13)

subplot(232)
plot(data(:,WaypointIndex), data(:,RoadHeightRMIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Road Height RM [m]')
set(gca,'FontSize',13)

subplot(233)
plot(data(:,WaypointIndex), data(:,RoadHeightRRIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Road Height RR [m]')
set(gca,'FontSize',13)

subplot(234)
plot(data(:,WaypointIndex), data(:,RoadHeightLFIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Road Height LF [m]')
set(gca,'FontSize',13)

subplot(235)
plot(data(:,WaypointIndex), data(:,RoadHeightLMIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Road Height LM [m]')
set(gca,'FontSize',13)

subplot(236)
plot(data(:,WaypointIndex), data(:,RoadHeightLRIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Road Height LR [m]')
set(gca,'FontSize',13)

figure
set(gcf,'Color',[1,1,1])
subplot(231)
plot(data(:,WaypointIndex), data(:,SuspensionRFIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Suspension Angle RF [m]')
set(gca,'FontSize',13)

subplot(232)
plot(data(:,WaypointIndex), data(:,SuspensionRMIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Suspension Angle RM [m]')
set(gca,'FontSize',13)

subplot(233)
plot(data(:,WaypointIndex), data(:,SuspensionRRIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Suspension Angle RR [m]')
set(gca,'FontSize',13)

subplot(234)
plot(data(:,WaypointIndex), data(:,SuspensionLFIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Suspension Angle LF [m]')
set(gca,'FontSize',13)

subplot(235)
plot(data(:,WaypointIndex), data(:,SuspensionLMIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Suspension Angle LM [m]')
set(gca,'FontSize',13)

subplot(236)
plot(data(:,WaypointIndex), data(:,SuspensionLRIndex),'LineWidth',2)
grid on
xlabel('WayPoint')
ylabel('Suspension Angle LR [m]')
set(gca,'FontSize',13)
% 
% figure('Name','RF')
% set(gcf,'Color',[1,1,1])
% subplot(311)
% plot(data(:,TimeIndex),data(:,RF_tire_x))
% grid on
% subplot(312)
% plot(data(:,TimeIndex),data(:,RF_tire_y))
% grid on
% subplot(313)
% plot(data(:,TimeIndex),data(:,RF_tire_z))
% grid on
% 
% figure('Name','RM')
% set(gcf,'Color',[1,1,1])
% subplot(311)
% plot(data(:,TimeIndex),data(:,RM_tire_x))
% grid on
% subplot(312)
% plot(data(:,TimeIndex),data(:,RM_tire_y))
% grid on
% subplot(313)
% plot(data(:,TimeIndex),data(:,RM_tire_z))
% grid on
% 
% figure('Name','RR')
% set(gcf,'Color',[1,1,1])
% subplot(311)
% plot(data(:,TimeIndex),data(:,RR_tire_x))
% grid on
% subplot(312)
% plot(data(:,TimeIndex),data(:,RR_tire_y))
% grid on
% subplot(313)
% plot(data(:,TimeIndex),data(:,RR_tire_z))
% grid on
% 
% figure('Name','LF')
% set(gcf,'Color',[1,1,1])
% subplot(311)
% plot(data(:,TimeIndex),data(:,LF_tire_x))
% grid on
% subplot(312)
% plot(data(:,TimeIndex),data(:,LF_tire_y))
% grid on
% subplot(313)
% plot(data(:,TimeIndex),data(:,LF_tire_z))
% grid on
% 
% figure('Name','LM')
% set(gcf,'Color',[1,1,1])
% subplot(311)
% plot(data(:,TimeIndex),data(:,LM_tire_x))
% grid on
% subplot(312)
% plot(data(:,TimeIndex),data(:,LM_tire_y))
% grid on
% subplot(313)
% plot(data(:,TimeIndex),data(:,LM_tire_z))
% grid on
% 
% figure('Name','LR')
% set(gcf,'Color',[1,1,1])
% subplot(311)
% plot(data(:,TimeIndex),data(:,LR_tire_x))
% grid on
% subplot(312)
% plot(data(:,TimeIndex),data(:,LR_tire_y))
% grid on
% subplot(313)
% plot(data(:,TimeIndex),data(:,LR_tire_z))
% grid on