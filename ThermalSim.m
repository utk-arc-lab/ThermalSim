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
			% length(kActivateTime)
			% kActivateTime

			% Configure Thermal Model
			thermal_model = ThermalSim.InitializeTransientThermalModel();
			% pde_mesh_nodes = ThermalSim.BuildSingleWallPDEMeshNodes(props);
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
			new_kActivateIndices = new_kActivateIndices + 1;
			kActivateTime(2:end) = kActivateTime(new_kActivateIndices);

			ThermalSim.SetInitialNodeTemperature(props,thermal_model,props.melt_temp,1:length(wall_faces));
			ThermalSim.SetFaceThermalConductivities(props,thermal_model,@NodeActivationFunction);

			% ThermalSim.ElementWiseSimSetup(thermal_sim_properties,thermal_model);

			% ThermalSim.SetAllButFirstInitialFaceTemperature(props,thermal_model);
			% ThermalSim.SetFaceThermalConductivities(props,thermal_model,@NodeActivationFunction);

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
			if(~isa(thermal_path,'ThermalPath'))
				fprintf('ThermalSim::BuildSingleWallPDEMeshNodesFromPath: Input not a ThermalPath\n');
				pde_mesh_nodes = [];
				return;
			end%if

			if(isempty(thermal_path.x) || isempty(thermal_path.y))
				fprintf('ThermalSim::BuildSingleWallPDEMeshNodesFromPath: ThermalPath not built!\n');
				pde_mesh_nodes = [];
				return;
			end%if
			
			node_length = thermal_path.wall_length / (thermal_path.nodes_per_layer - 1);
			node_height = thermal_path.bead_height;

			layers = zeros(10,length(thermal_path.x));

			for i = 1:length(thermal_path.x)
				x_i = thermal_path.x(i) + thermal_path.base_origin(1) + thermal_path.wall_origin(1);
				y_i = thermal_path.y(i) + thermal_path.base_origin(2) + thermal_path.wall_origin(2);

				x = [x_i - node_length / 2; x_i + node_length / 2; x_i + node_length / 2; x_i - node_length / 2];
				y = [y_i; y_i; y_i + node_height; y_i + node_height];

				layers(:,i) = ThermalSim.DecsgRectangle(x,y);
			end%for i

			base_x = [thermal_path.base_origin(1); thermal_path.base_origin(1) + thermal_path.base_length; ...
			 thermal_path.base_origin(1) + thermal_path.base_length; thermal_path.base_origin(1)];

			base_y = [thermal_path.base_origin(2); thermal_path.base_origin(2); ...
			 thermal_path.base_origin(2) + thermal_path.base_thick; thermal_path.base_origin(2) + thermal_path.base_thick;];

			base = ThermalSim.DecsgRectangle(base_x,base_y);

			pde_mesh_nodes = [base layers];

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
			% wall_elements = ThermalSim.GetSortedMeshElementsInRange(thermal_model,wall_x_range,wall_y_range);
		end%func GetSortedBuildElements

		function wall_elements = GetWallMeshElements(thermal_model,thermal_path)
			node_length = thermal_path.wall_length / (thermal_path.nodes_per_layer - 1);
			node_height = thermal_path.bead_height;

			wall_elements = cell(1);

			for i = 1:length(thermal_path.x)
				x_i = thermal_path.x(i) + thermal_path.base_origin(1) + thermal_path.wall_origin(1);
				y_i = thermal_path.y(i) + thermal_path.base_origin(2) + thermal_path.wall_origin(2);

				x = [x_i - node_length / 2, x_i + node_length / 2];
				y = [y_i, y_i + node_height];

				wall_elements{i} = thermal_model.Mesh.findElements('box',x,y);
			end%for i
		end%func GetSortedMeshElements

		function ElementWiseSimSetup(thermal_sim_properties,thermal_model)
			props = thermal_sim_properties; % shorthand

			ThermalSim.SetElementInitialTemperature(thermal_model,props.ambient_T,props.base_elements);
			ThermalSim.SetElementInitialTemperature(thermal_model,props.melt_temp,wall_elements);

			ThermalSim.SetElementThermalConductivities(thermal_sim_properties,thermal_model,props.k,base_elements);
			ThermalSim.SetElementThermalConductivities(thermal_sim_properties,thermal_model, ...
				@ActivationFunctions.NodeActivationFunction,wall_elements);
		end%func ElementWiseSimSetup

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

		function SetAllButFirstInitialFaceTemperature(thermal_sim_properties,thermal_model)
			props = thermal_sim_properties; % shorthand
			n_faces = thermal_model.Geometry.NumFaces; % shorthand

			fprintf('Setting all wall faces to melt temperature... ');
			tic;
			% Pre allocate all Faces other than the baseplate to melting temperature
			for i = 2:n_faces
			    thermalIC(thermal_model,props.melt_temp,'Face',i);
			end
			fprintf('%1.3fs\n',toc);
		end%func SetAllButFirstInitialFaceTemperature

		function kActivateTime = CalculateNodeActivationTimes(thermal_sim_properties)
			props = thermal_sim_properties; % shorthand
			t = props.thermal_path; % shorthand

			fprintf('Calculating Node Activation Times... ');
			tic;
			% Parameter Calculation
			kActivateTime = zeros(1, t.n_layers.*t.nodes_per_layer);
			% Skip Face 1 because it is baseplate and initial condition of room temprops...
			n = 2; 
			%initialize time for 1 second to give a 0 starting condition
			time = 1;

			for i = 0:t.n_layers-1
			    t_offset = t.layer_wait;
			    
			    for j = 1:t.nodes_per_layer
			        timepernode = (props.node_length./t.travel_speed);
			        if j == 1
			            time = t_offset+time;
			        else
			            time = time+timepernode;
			        end
			        kActivateTime(n) = time;
			        n = n+1;
			    end
			end
			kActivateTime(2:end) = kActivateTime(2:end) - t.layer_wait;

			fprintf('%1.3fs\n',toc);
		end%func CalculateNodeActivationTimes

		function SetFaceThermalConductivities(thermal_sim_properties,thermal_model,specific_heat_expression)
			props = thermal_sim_properties; % shorthand
			n_faces = thermal_model.Geometry.NumFaces;

			fprintf('Setting Node Properties... ');
			tic;

			thermalProperties(thermal_model,'ThermalConductivity', props.k,'MassDensity', props.density, ...
			                                'SpecificHeat', props.Cp);
			thermalIC(thermal_model,props.ambient_T,'Face',1);

			for i = 2:n_faces
			    thermalProperties(thermal_model,'Face',i,'ThermalConductivity',specific_heat_expression, ...
			        'MassDensity', props.density, 'SpecificHeat', props.Cp);
			end%for i

			fprintf('%1.3fs\n',toc);

		end%func SetFaceThermalConductivities

		function decsg_vector = DecsgRectangle(x,y)
			decsg_vector = [3; 4; x(1); x(2); x(3); x(4); y(1); y(2); y(3); y(4)];
		end%func DecsgRectangle

	end% Private Static Methods

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