clc, close all, clear all

% Set up custom properties
t = ThermalPath;
t.n_layers = 1;
t.wall_origin = [(t.base_length - t.wall_length)./2,t.base_thick];
t.nodes_per_layer = 50;

ThermalPathBuilder.BuildSolidWall(t);
% ThermalPathBuilder.BuildCastle(t);
% t.contours{end}.waypoints = flip(t.contours{end}.waypoints);

p = ThermalSimProperties(t); % default
% Options
p.bool_activate_nodes_iteratively = true;
p.bool_enable_internal_heat_generation = true;
% Parameters
p.weld_time_offset = 0.1; % s
p.baseplate_convection_coefficient = 1;

thermalmodel = ThermalSim.ConfigureSingleWallSim(p);

sim_time = 60; % s
sim_fps = 30;

[R,tlist] = ThermalSim.SolveSingleWallSim(thermalmodel,sim_time,sim_fps);
