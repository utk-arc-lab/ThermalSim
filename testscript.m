clc, close all, clear all

% Set up custom properties
t = ThermalPath;
t.n_layers = 3;
t.BuildPath();

p = ThermalSimProperties(t); % default

thermalmodel = ThermalSim.ConfigureSingleWallSim(p);
% [R,tlist] = ThermalSim.SolveSingleWallSim(thermalmodel,240,30);

% pdeplot(thermalmodel,'ElementLabels','on');

% load 4stest;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);

% SingleWallScript;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);