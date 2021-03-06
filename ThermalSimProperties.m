classdef ThermalSimProperties < handle
	properties % DEFAULT
		% Arg in
		thermal_path;

		% Bool config
		bool_activate_nodes_iteratively = true;
		bool_enable_internal_heat_generation = true; 

		% Set
		melt_temp = 1873.15; %Kelvin, based from above mild steel melt temps
		density = 7.87e-6; %low-carbon steel 7.87 g/cm3 converted to kg/mm3
		Cp = 620; % mild steel 620 J/kg K (Specific Heat)
		k = 49.6; % mild steel 1%C average vale W/m-K (Thermal Conductivity)
		ambient_T = 293.15; %Kelvin (20Deg C)
		Q_dot = 1200;
		
		weld_time_offset = 0.1; % s
		baseplate_convection_coefficient = 5.75*10^-6;

		% Calculated
		base_mass;
		node_width;
		node_length;
		node_thick;
		lump_mass;

		% Assigned in ThermalSim() for use later
		base_elements = [];
		wall_elements = [];
		base_faces = [];
		wall_faces = [];
		kActivateTime = [];
		kStartTime = [];
		kEndTime = [];
		
	end%properties

	methods
		function obj = ThermalSimProperties(thermal_path)
			if(isempty(thermal_path.contours))
				fprintf('ThermalSimProperties::ThermalSimProperties: Thermal Path not built! Use ThermalPathBuilder first!\n');
				return;
			end%if
			p = thermal_path; % shorthand
			
			% Set
			obj.thermal_path = p;

			% Calculate any calculated properties
			obj.base_mass = (p.base_width.*p.base_length.*p.base_thick).*obj.density; % kg
			obj.node_width = p.bead_width;
			obj.node_length = p.wall_length./p.nodes_per_layer; %% mm based on sim params
			obj.node_thick = p.bead_height;
			obj.lump_mass = (obj.node_width.*obj.node_length.*obj.node_thick).*obj.density; % kg
		end%func
	end%methods
end%func 