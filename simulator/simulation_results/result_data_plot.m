clc; clear all; close all;

data_indx = 1;
[data, txt] = xlsread(sprintf('%d/SimulationResult.csv', data_indx));

TimeIndex = 1;
VelocityIndex = 29;
WaypointIndex = 30;
DesVelocityIndex = 31;
SuspensionRFIndex = 14;
SuspensionRMIndex = 15;
SuspensionRRIndex = 16;
SuspensionLFIndex = 17;
SuspensionLMIndex = 18;
SuspensionLRIndex = 19;
Vd_ctrl_index = 31;
Vx_ctrl_index = 32;
RSM_index = 33;
PSM_index = 34;
LSM_index = 35;
VSM_index = 36;

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

figure
set(gcf,'Color',[1,1,1])
subplot(6,1,[1,2])
plot(data(:,TimeIndex), data(:,Vd_ctrl_index), 'LineWidth', 2)
hold on
plot(data(:,TimeIndex), data(:,Vx_ctrl_index), '--', 'LineWidth',2)
grid on
ylabel('Velocity [m/s]')
legend('Command','Velocity')
set(gca,'FontSize',13)

subplot(6,1,3)
plot(data(:,TimeIndex), data(:,RSM_index), 'LineWidth', 2)
grid on
ylabel('RSM')
ylim([0,1])
set(gca,'FontSize',13)

subplot(6,1,4)
plot(data(:,TimeIndex), data(:,PSM_index), 'LineWidth', 2)
grid on
ylabel('PSM')
ylim([0,1])
set(gca,'FontSize',13)

subplot(6,1,5)
plot(data(:,TimeIndex), data(:,LSM_index), 'LineWidth', 2)
grid on
ylabel('LSM')
ylim([0,1])
set(gca,'FontSize',13)

subplot(6,1,6)
plot(data(:,TimeIndex), data(:,VSM_index), 'LineWidth', 2)
grid on
ylabel('VSM')
ylim([0,1])
set(gca,'FontSize',13)
