clc; clear all; close all;

for i = 1 : 30
    [data, txt] = xlsread(sprintf('%d/SimulationResult.csv', i));
    vel(i,1) = mean(data(:,32));
    vel(i,2) = mean(data(:,32))*3600/1000;
end

mean(vel)