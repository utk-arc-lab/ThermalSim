clc, clear all, close all;
q_dot_list = [600,1200,6000,12000];

for q_dot = q_dot_list
	close all;

	% Set up custom properties
	t = ThermalPath;
	t.n_layers = 4;
	t.wall_origin = [(t.base_length - t.wall_length)./2,t.base_thick];
	t.nodes_per_layer = 20;

	ThermalPathBuilder.BuildSolidWall(t);
	% ThermalPathBuilder.BuildCastle(t);
	% t.contours{end}.waypoints = flip(t.contours{end}.waypoints);

	p = ThermalSimProperties(t); % default
	% Options
	p.bool_activate_nodes_iteratively = true;
	p.bool_enable_internal_heat_generation = true;
	% Parameters
	p.weld_time_offset = 0.1; % s
	p.Q_dot = q_dot;

	thermalmodel = ThermalSim.ConfigureSingleWallSim(p);

	sim_time = 360; % s
	sim_fps = 20;

	[R,tlist] = ThermalSim.SolveSingleWallSim(thermalmodel,sim_time,sim_fps);

	data_title = sprintf('ThermalSimQdot%1.3f',q_dot);
	ThermalPlotter.PDEAnimate(thermalmodel,R,tlist,true,data_title);

	save(data_title);
end%for q_dot