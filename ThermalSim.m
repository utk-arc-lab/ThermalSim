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

			p = thermal_sim_properties; % shorthand

			% Setup
			global kActivateTime;
			kActivateTime = ThermalSim.CalculateNodeActivationTimes(p);

			% Configure Thermal Model
			thermal_model = ThermalSim.InitializeTransientThermalModel();
			pde_mesh_nodes = ThermalSim.BuildSingleWallPDEMeshNodes(p);
			pde_mesh = ThermalSim.BuildSingleWallPDEMeshFromNodes(thermal_model,pde_mesh_nodes);

			[base_elements,wall_elements] = ThermalSim.GetSortedBuildElements(p,thermal_model);

			ThermalSim.SetAllButFirstInitialFaceTemperature(p,thermal_model);
			ThermalSim.SetFaceThermalConductivities(p,thermal_model,@NodeActivationFunction);

		end%func SingleWallSim
	end% Public Static Methods

	methods(Static, Access = 'private')
		function thermal_model = InitializeTransientThermalModel()
			fprintf('Building Thermal Model... ');
			tic;
			thermal_model = createpde('thermal','transient');
			fprintf('%1.3fs\n',toc);
		end%func InitializeTransientThermalModel

		function pde_mesh_nodes = BuildSingleWallPDEMeshNodes(thermal_sim_properties)
			p = thermal_sim_properties; % shorthand

			fprintf('Building PDE Mesh Nodes... ');
			tic;
			Base = [3 4 p.base_origin(1) p.base_length p.base_length p.base_origin(1) p.base_origin(2) p.base_origin(2) p.base_thick p.base_thick]';

			%Preallocate Array Size for Layers
			layers = zeros(10, p.num_layer.*p.nodeperlayer);
			n = 1;

			%Index Array Size Based on Node Numbers
			for i = 1:p.num_layer
			    bot_layerheight = ((i-1)*(p.node_thick))+p.wall_origin(2);
			    x_offset = p.wall_origin(1);
			    
			    for j = 0:p.nodeperlayer-1
			        layers(:,n) = [  3
			                                4
			                                x_offset+(j.*p.node_length)
			                                x_offset+((1+j).*p.node_length)
			                                x_offset+((1+j).*p.node_length)
			                                x_offset+(j.*p.node_length)
			                                bot_layerheight
			                                bot_layerheight
			                                bot_layerheight+p.node_thick
			                                bot_layerheight+p.node_thick];
			        
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
			p = thermal_sim_properties; % shorthand

			base_x_range = [0,p.base_length] + p.base_origin(1);
			base_y_range = [0,p.base_thick] + p.base_origin(2);
			wall_x_range = [0,(p.nodeperlayer * p.node_length)] + p.wall_origin(1);
			wall_y_range = [0,(p.num_layer * p.node_thick)] + p.wall_origin(2);

			% Query
			base_elements = ThermalSim.GetSortedMeshElementsInRange(thermal_model,base_x_range,base_y_range);
			wall_elements = ThermalSim.GetSortedMeshElementsInRange(thermal_model,wall_x_range,wall_y_range);
		end%func GetSortedBuildElements

		function sorted_elements = GetSortedMeshElementsInRange(thermal_model,x_range,y_range)
			% Extract
			mesh = thermal_model.Mesh;
			nodes = mesh.Nodes;
			elements = mesh.Elements;

			% Subset
			element_subset = mesh.findElements('box',x_range,y_range);
			
			% Preallocate
			[~,n_elements] = size(element_subset);
			element_centroids = zeros(2,n_elements);

			% Calculate Centroids
			for i = 1:n_elements
				element_nodes = elements(:,element_subset(i));

				nodes_x = nodes(1,element_nodes);
				nodes_y = nodes(2,element_nodes);

				nodes_avg_x = (max(nodes_x) + min(nodes_x)) / 2;
				nodes_avg_y = (max(nodes_y) + min(nodes_y)) / 2;
				element_centroids(:,i) = [nodes_avg_x;nodes_avg_y];
			end%for i

			[~,x_sort_indices] = sort(element_centroids(1,:));
			[~,y_sort_indices] = sort(element_centroids(2,:));



		end%func GetSortedMeshElements

		function SetAllButFirstInitialFaceTemperature(thermal_sim_properties,thermal_model)
			p = thermal_sim_properties; % shorthand
			n_faces = thermal_model.Geometry.NumFaces;

			fprintf('Setting all wall faces to melt temperature... ');
			tic;
			% Pre allocate all Faces other than the baseplate to melting temperature
			for i = 2:n_faces
			    thermalIC(thermal_model,p.melt_temp,'Face',i);
			end
			fprintf('%1.3fs\n',toc);
		end%func SetAllButFirstInitialFaceTemperature

		function kActivateTime = CalculateNodeActivationTimes(thermal_sim_properties)
			p = thermal_sim_properties; % shorthand

			fprintf('Calculating Node Activation Times... ');
			tic;
			% Parameter Calculation
			kActivateTime = zeros(1, p.num_layer.*p.nodeperlayer);
			% Skip Face 1 because it is baseplate and initial condition of room temp...
			n = 2; 
			%initialize time for 1 second to give a 0 starting condition
			time = 1;

			for i = 0:p.num_layer-1
			    t_offset = p.layerwait;
			    
			    for j = 1:p.nodeperlayer
			        timepernode = (p.node_length./p.travel_speed);
			        if j == 1
			            time = t_offset+time;
			        else
			            time = time+timepernode;
			        end
			        kActivateTime(n) = time;
			        n = n+1;
			    end
			end
			kActivateTime(2:end) = kActivateTime(2:end) - p.layerwait;

			fprintf('%1.3fs\n',toc);
		end%func CalculateNodeActivationTimes

		function SetFaceThermalConductivities(thermal_sim_properties,thermal_model,specific_heat_expression)
			p = thermal_sim_properties; % shorthand
			n_faces = thermal_model.Geometry.NumFaces;

			fprintf('Setting Node Properties... ');
			tic;

			thermalProperties(thermal_model,'ThermalConductivity', p.k,'MassDensity', p.density, ...
			                                'SpecificHeat', p.Cp);
			thermalIC(thermal_model,p.ambient_T,'Face',1);

			for i = 2:n_faces
			    thermalProperties(thermal_model,'Face',i,'ThermalConductivity',specific_heat_expression, ...
			        'MassDensity', p.density, 'SpecificHeat', p.Cp);
			end%for i

			fprintf('%1.3fs\n',toc);

		end%func SetFaceThermalConductivities
	end% Private Static Methods

end%class ThermalSim