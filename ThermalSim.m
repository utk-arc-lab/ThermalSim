classdef ThermalSim
	methods(Static)
		function [results,tlist] = SolveSingleWallSim(thermal_model,sim_time,sample_fps)
			fprintf('Running Simulation... ');
			tic;

			tlist = linspace(0,sim_time,sim_time*sample_fps);

			results = solve(thermal_model,tlist);

			fprintf('%1.3fs\n',toc);
			return;
		end%func SolveSingleWallSim

		function thermal_model = ConfigureSingleWallSim(thermal_sim_properties)
			if(~isa(thermal_sim_properties,'ThermalSimProperties'))
				fprintf('ThermalSim::SingleWallSim: Input 1 not a ThermalSimProperties!\n');
				results = [];
				return;
			end%if

			props = thermal_sim_properties; % shorthand

			% Setup
			global kActivateTime;
			kActivateTime = ThermalSim.CalculateNodeActivationTimes(props);

			% Configure Thermal Model
			thermal_model = ThermalSim.InitializeTransientThermalModel();
			pde_mesh_nodes = ThermalSim.BuildSingleWallPDEMeshNodesFromPath(props.thermal_path);
			pde_mesh = ThermalSim.BuildSingleWallPDEMeshFromNodes(thermal_model,pde_mesh_nodes);

			[base_elements,wall_elements] = ThermalSim.GetSortedBuildElements(props,thermal_model);

			% Set
			props.base_elements = base_elements;
			props.wall_elements = wall_elements;

			[base_faces,wall_faces] = ThermalSim.GetBaseAndWallFaces(props,thermal_model);

			props.base_faces = base_faces;
			props.wall_faces = wall_faces;

			[~,new_kActivateIndices] = sort(wall_faces);
			new_kActivateIndices = new_kActivateIndices;
			kActivateTime = kActivateTime(new_kActivateIndices);

			global kStartTime;
			kStartTime = kActivateTime + 0.1;

			global kEndTime;
			kEndTime = ThermalSim.CalculateNodeEndTimes(props,kActivateTime);

			base_activation_time = 0;
			kActivateTime = [kActivateTime(1:base_faces - 1),base_activation_time,kActivateTime(base_faces:end)];
			kEndTime = [kEndTime(1:base_faces - 1), base_activation_time, kEndTime(base_faces:end)];
			kStartTime = [kStartTime(1:base_faces - 1), base_activation_time, kStartTime(base_faces:end)];

			ThermalSim.SetInitialNodeTemperature(props,thermal_model,props.ambient_T,1:length(wall_faces));
			% ThermalSim.SetFaceThermalConductivities(props,thermal_model,@ActivationFunctions.NodeActivationFunction);
			ThermalSim.SetFaceThermalConductivities(props,thermal_model,props.k);

			ThermalSim.SetFaceInternalHeatGeneration(props,thermal_model,@ActivationFunctions.InternalHeatingFunction);

			thermalBC(thermal_model,'Edge',1,'ConvectionCoefficient',5.75,'AmbientTemperature',props.ambient_T);
		end%func SingleWallSim
	end% Public Static Methods

	methods(Static, Access = 'private')
		function thermal_model = InitializeTransientThermalModel()
			fprintf('Building Thermal Model... ');
			tic;
			thermal_model = createpde('thermal','transient');
			fprintf('%1.3fs\n',toc);
		end%func InitializeTransientThermalModel

		function  pde_mesh_nodes = BuildSingleWallPDEMeshNodesFromPath(thermal_path)
			node_length = thermal_path.node_length;
			node_height = thermal_path.bead_height;

			layers_decsg = zeros(10,1);
			k = 0;

			for i = 1:length(thermal_path.contours)
				[x,y] = thermal_path.contours{i}.PathVectors();

				for j = 1:length(x)
					k = k + 1;
					x_i = x(j) + thermal_path.base_origin(1) + thermal_path.wall_origin(1);
					y_i = y(j) + thermal_path.base_origin(2) + thermal_path.wall_origin(2);

					x_node = [x_i - node_length / 2; x_i + node_length / 2; x_i + node_length / 2; x_i - node_length / 2];
					y_node = [y_i; y_i; y_i + node_height; y_i + node_height];

					layers_decsg(:,k) = ThermalSim.DecsgRectangle(x_node,y_node);
				end%for j
			end%for i

			base_x = [thermal_path.base_origin(1); thermal_path.base_origin(1) + thermal_path.base_length; ...
			 thermal_path.base_origin(1) + thermal_path.base_length; thermal_path.base_origin(1)];

			base_y = [thermal_path.base_origin(2); thermal_path.base_origin(2); ...
			 thermal_path.base_origin(2) + thermal_path.base_thick; thermal_path.base_origin(2) + thermal_path.base_thick;];

			base_decsg = ThermalSim.DecsgRectangle(base_x,base_y);

			pde_mesh_nodes = [base_decsg layers_decsg];

		end%func BuildSingleWallPDEMeshNodesFromPath

		function pde_mesh_nodes = BuildSingleWallPDEMeshNodes(thermal_sim_properties)
			props = thermal_sim_properties; % shorthand
			t = props.thermal_path; % shorthand

			fprintf('Building PDE Mesh Nodes... ');
			tic;
			Base = [3 4 t.base_origin(1) t.base_length t.base_length t.base_origin(1) t.base_origin(2) t.base_origin(2) t.base_thick t.base_thick]';

			%Preallocate Array Size for Layers
			layers = zeros(10, t.n_layers.*t.nodes_per_layer);
			n = 1;

			%Index Array Size Based on Node Numbers
			for i = 1:t.n_layers
			    bot_layerheight = ((i-1)*(props.node_thick))+t.wall_origin(2);
			    x_offset = t.wall_origin(1);
			    
			    for j = 0:t.nodes_per_layer-1
			        layers(:,n) = [  3
			                                4
			                                x_offset+(j.*props.node_length)
			                                x_offset+((1+j).*props.node_length)
			                                x_offset+((1+j).*props.node_length)
			                                x_offset+(j.*props.node_length)
			                                bot_layerheight
			                                bot_layerheight
			                                bot_layerheight+props.node_thick
			                                bot_layerheight+props.node_thick];
			        
			        n = n+1;
			    end%for j

			end%for i

			pde_mesh_nodes = [Base layers];

			fprintf('%1.3fs\n',toc);
		end%func BuildSingleWallPDEMeshNodes

		function pde_mesh = BuildSingleWallPDEMeshFromNodes(thermal_model,pde_mesh_nodes)
			fprintf('Generating PDE Mesh... ');
			dl = decsg(pde_mesh_nodes);
			pde_mesh = geometryFromEdges(thermal_model,dl);
			generateMesh(thermal_model,'Hmin',20);
			fprintf('%1.3fs\n',toc);
		end%func BuildSingleWallPDEMeshFromNodes

		function [base_elements,wall_elements] = GetSortedBuildElements(thermal_sim_properties,thermal_model)
			props = thermal_sim_properties; % shorthand
			t = props.thermal_path; % shorthand

			base_x_range = [0,t.base_length] + t.base_origin(1);
			base_y_range = [0,t.base_thick] + t.base_origin(2);
			wall_x_range = [0,(t.nodes_per_layer * props.node_length)] + t.wall_origin(1);
			wall_y_range = [0,(t.n_layers * props.node_thick)] + t.wall_origin(2);

			% Query
			base_elements = thermal_model.Mesh.findElements('box',base_x_range,base_y_range);
			wall_elements = ThermalSim.GetWallMeshElements(thermal_model,t);
		end%func GetSortedBuildElements

		function wall_elements = GetWallMeshElements(thermal_model,thermal_path)
			node_length = thermal_path.node_length;
			node_height = thermal_path.bead_height;

			wall_elements = cell(1);
			k = 0;

			for i = 1:length(thermal_path.contours)
				[x,y] = thermal_path.contours{i}.PathVectors();

				for j = 1:length(x)
					k = k + 1;

					x_i = x(j) + thermal_path.base_origin(1) + thermal_path.wall_origin(1);
					y_i = y(j) + thermal_path.base_origin(2) + thermal_path.wall_origin(2);

					x_range = [x_i - node_length / 2, x_i + node_length / 2];
					y_range = [y_i, y_i + node_height];

					wall_elements{k} = thermal_model.Mesh.findElements('box',x_range,y_range);
				end%for j
			end%for i
		end%func GetSortedMeshElements

		function SetElementInitialTemperature(thermal_model,temperature,elements)
			for i = 1:length(elements)
				thermalIC(thermal_model,temperature,'Element',i);
			end%for i
		end%func SetElementInitialTemperature

		function SetElementThermalConductivities(thermal_sim_properties,thermal_model,specific_heat_expression,elements)
			props = thermal_sim_properties; % shorthand			

			for i = 1:length(elements)
				thermalProperties(thermal_model,'ThermalConductivity',specific_heat_expression, ...
					'MassDensity', props.density, 'SpecificHeat', props.Cp)
			end%for i
		end%func SetElementThermalConductivities

		function SetInitialNodeTemperature(thermal_sim_properties,thermal_model,temperature,node_indices)
			props = thermal_sim_properties; % shorthand
			faces = props.wall_faces;

			for current_node = node_indices
				thermalIC(thermal_model,temperature,'face',faces(current_node));
			end%for node_indices
		end%func SetInitialNodeTemperatures

		function kActivateTime = CalculateNodeActivationTimes(thermal_sim_properties)
			props = thermal_sim_properties; % shorthand
			t = props.thermal_path; % shorthand

			fprintf('Calculating Node Activation Times... ');
			tic;

			% Parameter Calculation
			kActivateTime = [];

			k = 0; 
			%initialize time for 1 second to give a 0 starting condition
			time = 1;

			for i = 1:length(t.contours)
			    current_contour = t.contours{i};
			    t_offset = current_contour.wait_after_time;
			    
			    for j = 1:length(current_contour.waypoints)
			    	k = k+1;
			        node_time = (props.node_length./current_contour.waypoints{i}.travel_speed);
			        
			        if j == 1
			            time = t_offset+time;
			        else
			            time = time+node_time;
			        end
			        kActivateTime(k) = time;
			    end%for j
			end%for i

			kActivateTime = kActivateTime - t.contours{1}.wait_after_time;

			fprintf('%1.3fs\n',toc);
		end%func CalculateNodeActivationTimes

		function kEndTime = CalculateNodeEndTimes(thermal_sim_properties, kActivateTime)
			props = thermal_sim_properties; % shorthand
			t = props.thermal_path;

			kEndTime = kActivateTime;
			k = 0;

			for i = 1:length(t.contours)
				current_contour = t.contours{i};
				
				for j = 1:length(current_contour.waypoints)
					k = k + 1;
					kEndTime(k) = kEndTime(k) + (props.node_length / current_contour.waypoints{i}.travel_speed);
				end%for j
			end%for i
		end%func CalculateNodeEndTimes

		function SetFaceThermalConductivities(thermal_sim_properties,thermal_model,specific_heat_expression)
			props = thermal_sim_properties; % shorthand
			n_faces = thermal_model.Geometry.NumFaces;

			fprintf('Setting Node Properties... ');
			tic;

			% Set all faces to the default properties
			thermalProperties(thermal_model,'ThermalConductivity', props.k,'MassDensity', props.density, ...
			                                'SpecificHeat', props.Cp);

			% Set base face to ambient T
			thermalIC(thermal_model,props.ambient_T,'Face',props.base_faces);

			% Update wall face properties to include specific heat expression
			for current_face = props.wall_faces
			    thermalProperties(thermal_model,'Face',current_face,'ThermalConductivity',specific_heat_expression, ...
			        'MassDensity', props.density, 'SpecificHeat', props.Cp);
			end%for i

			fprintf('%1.3fs\n',toc);

		end%func SetFaceThermalConductivities

		function SetFaceInternalHeatGeneration(thermal_sim_properties,thermal_model,internal_heat_generation_expression)
			props = thermal_sim_properties; % shorthand

			for current_face = props.wall_faces
				internalHeatSource(thermal_model,internal_heat_generation_expression,'Face',current_face);
			end%for current face
		end%func SetFaceInternalHeatGeneration

		function decsg_vector = DecsgRectangle(x,y)
			decsg_vector = [3; 4; x(1); x(2); x(3); x(4); y(1); y(2); y(3); y(4)];
		end%func DecsgRectangle

	end% Private Static Methods

	% ======================= %
	% TESTING FUNCTIONS BELOW %
	% ======================= %

	methods(Static)
		function [base_face,wall_faces] = GetBaseAndWallFaces(thermal_sim_properties,thermal_model)
			props = thermal_sim_properties; % shorthand
			geometry = thermal_model.Geometry;

			wall_faces = [];
			for i = 1:thermal_model.Geometry.NumFaces
				[x_bounds,y_bounds] = ThermalSim.GetFaceBounds(thermal_model,i);
				elements_in_face = thermal_model.Mesh.findElements('box',x_bounds,y_bounds);
				elements_in_face = sort(elements_in_face);

				bool_node_found = false;
				for j = 1:length(props.wall_elements)
					current_elements = sort(props.wall_elements{j});

					if(size(current_elements) == size(elements_in_face))

						if(current_elements == elements_in_face)
							bool_node_found = true;
							wall_faces(j) = i;
							break;
						end%if

					end%if
				end%for j
				
				if(~bool_node_found)
					base_face = i;
				end%if
			end%for i
		end%func GetBaseAndWallFaces

		function [x_bounds,y_bounds] = GetElementBounds(thermal_model,elements)
			x = [];
			y = x;

			for current_element = elements
				node_set = thermal_model.Mesh.Elements(:,current_element)';
				
				for i = 1:length(node_set)
					current_node = node_set(i);
					x = [x thermal_model.Mesh.Nodes(1,current_node)];
					y = [y thermal_model.Mesh.Nodes(2,current_node)];
				end%for i

				x_bounds = [min(x),max(x)];
				y_bounds = [min(y),max(y)];
			end%for elements
		end%func GetElementBounds

		function [x_bounds,y_bounds] = GetFaceBounds(thermal_model,face_index)
			geometry = thermal_model.Geometry.geom;

			face_logical = geometry(7,:) == face_index;

			x = geometry(2:3,face_logical);
			y = geometry(4:5,face_logical);

			x_bounds = [min(x(:)),max(x(:))];
			y_bounds = [min(y(:)),max(y(:))];
		end%func GetFaceBounds
	end%static methods TEST

end%class ThermalSim