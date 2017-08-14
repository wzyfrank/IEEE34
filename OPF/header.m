clear all;
clc;
close all;

format short g;
% import Zbus and Cbus
Zbus = load('Zbus.mat');
Zbus = Zbus.Zbus;
Cbus = load('Cbus.mat');
Cbus = Cbus.Cbus;

% init the voltage profile
[num,txt,raw] = xlsread('loads data.csv');
N_loads = size(num, 1);

Voltage_profile = ones(N_loads, 6);
Voltage_profile(:, 2) = 0;
Voltage_profile(:, 4) = -120;
Voltage_profile(:, 6) = 120;

%% 1. build a phase->index hash map
[num_bus,txt_bus,raw_bus] = xlsread('ieee34_EXP_Y.CSV');
phaseMap = containers.Map;
for k = 1:size(raw_bus, 1)
    if(isnumeric(raw_bus{k,1}))
        phaseMap(num2str(raw_bus{k,1})) = k;
    else
        phaseMap(raw_bus{k,1}) = k;
    end
end

Vbase = 24900 / sqrt(3);
% start iteration
for iter = 1:20
    % output loads
    % loads = load_process(Voltage_profile, phaseMap);
    loads = load_process(Voltage_profile, phaseMap, Vbase);
    [Voltage_profile, phasors] = ieee34_iter2(loads, Zbus, Cbus);
    % [Voltage_profile, phasors] = ieee34_hybrid2(loads, Zbus, Cbus);
    % update the loads
end

