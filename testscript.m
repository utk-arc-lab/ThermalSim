clc, close all, clear all

% Set up custom properties
p = ThermalSimProperties; % default
p.num_layer = 3;

thermalmodel = ThermalSim.ConfigureSingleWallSim(p);
% [R,tlist] = ThermalSim.SolveSingleWallSim(thermalmodel,240,30);

% pdeplot(thermalmodel,'ElementLabels','on');

% load 4stest;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);

% SingleWallScript;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);