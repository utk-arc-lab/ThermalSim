classdef Waypoint
	properties
		x;
		y;
		travel_speed;
	end%properties
	methods
		function obj = Waypoint(x,y,travel_speed)
			obj.x = x;
			obj.y = y;
			obj.travel_speed = travel_speed;
		end%Constructor
	end%methods
end%class Waypoint