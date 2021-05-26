classdef ThermalPathBuilder
	properties(Constant)
		default_travel_speed = 2.75; % mm/s
		default_layer_wait = 30; % s
	end%const

	methods(Static)
		function WEEEEE(thermal_path)

		end%func

		function BuildSolidWall(thermal_path)
			thermal_path.contours = ThermalPathBuilder.GenerateSolidWallContours(thermal_path);
		end%func BuildSolidWall

		function BuildCastle(thermal_path)
			contours = ThermalPathBuilder.GenerateSolidWallContours(thermal_path);
			
			top_layer = contours{end};
			n_points = length(top_layer.waypoints);

			castle_flags = [ceil(n_points / 4), floor(3*n_points / 4)];
			castle_waypoint_indices = [1:castle_flags(1), castle_flags(2):n_points];

			new_waypoints = cell(1);

			for i = 1:length(castle_waypoint_indices)
				old_point = top_layer.waypoints{castle_waypoint_indices(i)};
				new_waypoints{i} = Waypoint(old_point.x,old_point.y,old_point.travel_speed);
			end%for i

			contours{end}.waypoints = new_waypoints;
			thermal_path.contours = contours;
		end%func BuildCastle
	end%static methods

	methods(Static, Access = 'private')
		function waypoints = LinspaceWaypoints(x_range,y_range,n_waypoints)
			waypoints = cell(1);

			x = linspace(x_range(1),x_range(2),n_waypoints);
			y = linspace(y_range(1),y_range(2),n_waypoints);

			for i = 1:n_waypoints
				waypoints{i} = Waypoint(x(i),y(i),ThermalPathBuilder.default_travel_speed);
			end%for i
		end%func LinspaceWaypoints

		function contours = GenerateSolidWallContours(thermal_path)
			layer_height = thermal_path.bead_height;
			n_layers = thermal_path.n_layers;
			wall_width = thermal_path.wall_length;
			nodes_per_layer = thermal_path.nodes_per_layer;
			node_length = thermal_path.node_length;

			contours = cell(1);

			for i = 1:n_layers
				current_layer_y = (i-1) * layer_height;

				x_range = [0,wall_width];
				y_range = [current_layer_y,current_layer_y];

				waypoints = ThermalPathBuilder.LinspaceWaypoints(x_range,y_range,nodes_per_layer);
				contours{i} = Contour(waypoints,ThermalPathBuilder.default_layer_wait);
			end%for i
		end%func GenerateSolidWallContours
	end%private methods
end%class ThermalPathBuilder