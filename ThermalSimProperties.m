classdef ThermalSimProperties
	properties % DEFAULT
		travel_speed = 2.75; %mm/sec
		melt_temp = 1873.15; %Kelvin, based from above mild steel melt temps

		bead_height = 2.25; %mm
		bead_width = 6; %mm

		density = 7.87e-6; %low-carbon steel 7.87 g/cm3 converted to kg/mm3
		Cp = 620; % mild steel 620 J/kg K (Specific Heat)
		k = 49.6; % mild steel 1%C average vale W/m-K (Thermal Conductivity)
		ambient_T = 293.15; %Kelvin (20Deg C)

		num_layer = 30; % 30 layer wall 
		nodeperlayer = 20; % Change this arbitrarily
		wall_length = 152.4; % 6 in long
		layerwait = 30; % seconds

		base_width = 152.4; %mm
		base_length = 304.8; %mm
		base_thick = 12.7; %mm

		base_origin = [0,0];

		% Calculated
		base_mass;
		node_width;
		node_length;
		node_thick;
		lump_mass;
		wall_origin;
	end%properties

	methods
		function obj = ThermalSimProperties
			% Calculate any calculated properties
			obj.base_mass = (obj.base_width.*obj.base_length.*obj.base_thick).*obj.density; % kg
			obj.node_width = obj.bead_width;
			obj.node_length = obj.wall_length./obj.nodeperlayer; %% mm based on sim params
			obj.node_thick = obj.bead_height;
			obj.lump_mass = (obj.node_width.*obj.node_length.*obj.node_thick).*obj.density; % kg
			obj.wall_origin = [(obj.base_length - obj.wall_length)./2,obj.base_thick];
		end%func
	end%methods
end%func 