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
            % thermal_conductivity = 0.00001;
            thermal_conductivity = 0;
        end%if
        
        thermal_conductivities(i) = thermal_conductivity;
    end%for i
end%func