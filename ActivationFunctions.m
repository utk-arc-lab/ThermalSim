classdef ActivationFunctions
	methods(Static)
		function thermal_conductivities = NodeActivationFunction(location,state)
		    global kActivateTime;

		    if any(isnan(state.time))
		        thermal_conductivities = nan(size(state.u));
		        return;
		    end
		    face_indices = location.subdomain;
		    
		    thermal_conductivities = zeros(size(face_indices));
		    for i = 1:length(thermal_conductivities)
		        current_face = face_indices(i);
		        
		        if(state.time > kActivateTime(current_face))
		            thermal_conductivity = 49.8;
		        else
		            thermal_conductivity = 0;
		        end%if
		        
		        thermal_conductivities(i) = thermal_conductivity;
		    end%for i
		end%func NodeActivationFunction

		function internal_heat_generation = InternalHeatingFunction(location,state)
			global kStartTime;
			global kEndTime;
			global Q_dot;

		    if any(isnan(state.time))
		        internal_heat_generation = nan(size(state.u));
		        return;
		    end%if
		    face_indices = location.subdomain;
		    
		    internal_heat_generation = zeros(size(face_indices));
		    for i = 1:length(internal_heat_generation)
		        current_face = face_indices(i);
		        
		        start_time = kStartTime(current_face);

		        end_time = kEndTime(current_face);
		        
		        if(state.time > start_time && state.time < end_time)
		            internal_heat = Q_dot;
		        else
		            % thermal_conductivity = 0.00001;
		            internal_heat = 0;
		        end%if
		        
		        internal_heat_generation(i) = internal_heat;
		    end%for i
		end%func InternalHeatingFunction
	end%Static methods
end%class