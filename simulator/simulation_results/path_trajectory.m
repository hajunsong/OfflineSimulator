clc; clear all; close all;

simulation_indx = 1;
addpath(sprintf('%d',simulation_indx));

data = csvread('SimulationResult.csv',1,0); % 시뮬레이션 결과 데이터

east = load('../east_temp.csv');
north = load('../north_temp.csv');

figure
plot(east,north,'LineWidth',2)
hold on
plot(data(:,2), data(:,3),'--','LineWidth',2)
grid on