classdef ThermalPlotter
	properties(Constant)
		plot_colormap = jet;
	end%const

	methods(Static)
		function PDEAnimate(thermal_model,results,time_vector,bool_render,file_name)
			if(bool_render)
				if(nargin == 5)
					video_filename = strcat('Renders\',file_name);
				else
					video_filename = sprintf('Renders\ThermalSim%is',ceil(time_vector(end)));
				end%if
				v = VideoWriter(video_filename,'MPEG-4');
				v.FrameRate = 60;
				open(v);
			end%if

			f = figure('units','normalized','position',[0.1,0.1,0.8,0.8]);
			a = axes;

			[n_triangles,n_time_steps] = size(results.Temperature);

			for i = 1:n_time_steps
				graphics_handle = ThermalPlotter.PDEPlot(thermal_model,results,time_vector,i);
				% ThermalPlotter.CustomPDEPlot(thermal_model,results,time_vector,i);
				
				if(bool_render)
					frame=getframe(f);
					writeVideo(v,frame);
				else
					pause(0.01);
				end%if
			end%for i

			if(bool_render)
				close(v);
			end%if

		end%func PDEAnimate

		function PDECustomAnimate(thermal_sim_properties,thermal_model,results,time_vector,bool_render,file_name)
			if(bool_render)
				if(nargin == 6)
					video_filename = strcat('Renders\',file_name);
				else
					video_filename = sprintf('Renders\ThermalSim%is',ceil(time_vector(end)));
					title(sprintf('Time: %1.2fs',time_vector(time_index)));
				end%if
				v = VideoWriter(video_filename,'MPEG-4');
				v.FrameRate = 60;
				open(v);
			end%if

			f = figure('units','normalized','position',[0.1,0.1,0.8,0.8]);

			[n_triangles,n_time_steps] = size(results.Temperature);

			for i = 1:n_time_steps
				if i == 1
					a = ThermalPlotter.CustomPDEPlot(thermal_sim_properties,thermal_model,results,time_vector,i);
				else
					a.CData = results.Temperature(:,i);
				end%if

				title(sprintf('Time: %1.3fs',time_vector(i)));
				
				if(bool_render)
					frame=getframe(f);
					writeVideo(v,frame);
				else
					pause(0.01);
				end%if
			end%for i

			if(bool_render)
				close(v);
			end%if

		end%func PDEAnimate

		function h = PDEPlot(thermal_model,results,time_vector,time_index)
			temp = results.Temperature;
			% t_max = max(temp(:));
			t_max = 2000; % should be max base plate temperature
			t_min = 293;%min(temp(:));

			h = pdeplot(thermal_model,'XYData',temp(:,time_index),'colormap',ThermalPlotter.plot_colormap);
			
			grid on;
			colorbar;
			caxis([t_min,t_max]);
			xlabel('X (mm)');
			ylabel('Y (mm)');
			title(sprintf('Time: %1.2fs',time_vector(time_index)));
			ylim([0,50]);
			xlim([0,300]);
		end%func

		function ref = CustomPDEPlot(thermal_sim_properties,thermal_model,results,time_vector,time_index)
			node_list = thermal_model.Mesh.Nodes;
			temp = results.Temperature(:,time_index);

			temp_active_index = temp > 0;
			node_sub_list = node_list(:,temp_active_index);
			temp_sub = temp(temp_active_index);

			temp_range = [293,max(results.Temperature(:))];
			[~,n_points] = size(node_list);

			x = node_sub_list(1,:);
			y = node_sub_list(2,:);

			ref = scatter(x,y,120.*ones(1,length(temp_sub)),temp_sub,'filled');
			colormap jet;
			colorbar;
			caxis(temp_range);

			grid on;
			xlabel('X (mm)');
			ylabel('Y (mm)');
			title(sprintf('Time: %1.3fs',time_vector(time_index)));
			ylim([0,50]);
			xlim([0,300]);
		end%func
	end%static methods
end%class