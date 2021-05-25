clc, close all, clear all

% Set up custom properties
t = ThermalPath;
t.n_layers = 2;
t.BuildPath();
t.x(21:end) = flip(t.x(21:end));
% t.x = [t.x(1:25) t.x(35:end)];
% t.y = [t.y(1:25) t.y(35:end)];

p = ThermalSimProperties(t); % default

thermalmodel = ThermalSim.ConfigureSingleWallSim(p);
% [R,tlist] = ThermalSim.SolveSingleWallSim(thermalmodel,120,30);

% pdeplot(thermalmodel,'ElementLabels','on');

% load 4stest;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);

% SingleWallScript;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);