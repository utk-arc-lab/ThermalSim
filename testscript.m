clc, close all, clear all

% Set up custom properties
t = ThermalPath;
t.n_layers = 2;
t.wall_origin = [(t.base_length - t.wall_length)./2,t.base_thick];

% ThermalPathBuilder.BuildSolidWall(t);
ThermalPathBuilder.BuildCastle(t);
t.contours{end}.waypoints = flip(t.contours{end}.waypoints);

p = ThermalSimProperties(t); % default

thermalmodel = ThermalSim.ConfigureSingleWallSim(p);
[R,tlist] = ThermalSim.SolveSingleWallSim(thermalmodel,60,30);






% pdeplot(thermalmodel,'ElementLabels','on');

% load 4stest;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);

% SingleWallScript;
% ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true);