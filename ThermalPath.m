classdef ThermalPath < handle
	properties
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

		% Wall Properties
		wall_origin = [0,0];
		
		% Path
		contours;

		% Calculated
		node_length;

	end%properties

	methods
		function obj = ThermalPath()
			fprintf('ThermalPath::ThermalPath: Empty Thermal Path Initialized.\n');
			SetNodeLength(obj);
		end%func Constructor

		function obj = set.wall_length(obj,wall_length)
			obj.wall_length = wall_length;
			SetNodeLength(obj);
		end%func

		function obj = set.nodes_per_layer(obj,nodes_per_layer)
			obj.nodes_per_layer = nodes_per_layer;
			SetNodeLength(obj);
		end%func

		function obj = SetNodeLength(obj)
			obj.node_length = obj.wall_length / (obj.nodes_per_layer - 1);
		end%func
	end%methods
end%class ThermalPathSpecification