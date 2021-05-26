classdef Contour < handle
	properties
		waypoints;
		wait_after_time = 30;
	end%properties
	methods
		function obj = Contour(waypoints,wait_after_time)
			obj.waypoints = waypoints;
			obj.wait_after_time = wait_after_time;
		end%Constructor

		function [x,y] = PathVectors(obj)
			x = zeros(1,length(obj.waypoints));
			y = x;

			for i = 1:length(obj.waypoints)
				current_waypoint = obj.waypoints{i};
				x(i) = current_waypoint.x;
				y(i) = current_waypoint.y;
			end%for i
		end%func PathVectors
	end%methods
end%class Contour