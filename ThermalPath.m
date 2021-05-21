classdef ThermalPath < handle
	properties
		% Arg in
		layer_wait = 30; % seconds
		travel_speed = 2.75; % mm/sec

		% Set
		% Bead Properties
		bead_height = 2.25; %mm
		bead_width = 6; %mm

		% Build Properties
		n_layers = 30; % 30 layer wall 
		nodes_per_layer = 20; % Change this arbitrarily
		wall_length = 152.4; % 6 in long

		% Base Properties
		base_width = 152.4; %mm
		base_length = 304.8; %mm
		base_thick = 12.7; %mm
		base_origin = [0,0];

		% Calculated
		wall_origin;
		x; % mm
		y; % mm

	end%properties

	methods
		function obj = ThermalPath()
			% BuildPath(obj)
		end%func Constructor

		function BuildPath(obj)
			obj.wall_origin = [(obj.base_length - obj.wall_length)./2,obj.base_thick];

			x = [];
			y = [];
			for i = 1:obj.n_layers
				x = [x, linspace(0,obj.wall_length,obj.nodes_per_layer)];
				y = [y, ((i-1)*obj.bead_height).*ones(1,obj.nodes_per_layer)];
			end%for i
			obj.x = x;
			obj.y = y;
		end%func 
	end%methods
end%class ThermalPathSpecification