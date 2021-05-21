%% New Thermal Example Wall Single Layer Example

close all; clear all; clc;
warning('off','MATLAB:subscripting:noSubscriptsSpecified');

%% Necessary Planning Values for Simulation
fprintf('Initializing Parameters... ');
tic;
TS = 2.75; %mm/sec
% TS = 30; % test
melt_temp = 1873.15; %Kelvin, based from above mild steel melt temps

bead_height = 2.25; %mm
bead_width = 6; %mm

density = 7.87e-6; %low-carbon steel 7.87 g/cm3 converted to kg/mm3
Cp = 620; %mild steel 620 J/kg K (Specific Heat)
k = 49.6; %mild steel 1%C average vale W/m-K (Thermal Conductivity)
ambient_T = 293.15; %Kelvin (20Deg C)

num_layer = 1; % 30 layer wall 
nodeperlayer = 20; % Change this arbitrarily
wall_length = 152.4; %6 in long
layerwait = 1; %seconds
% layerwait = 0; % test

base_width = 152.4; %mm
base_length = 304.8; %mm
base_thick = 12.7; %mm
base_mass = (base_width.*base_length.*base_thick).*density; % kg

node_width = bead_width;
node_length = wall_length./nodeperlayer; %% mm based on sim params
node_thick = bead_height;
lump_mass = (node_width.*node_length.*node_thick).*density; % kg
fprintf('%1.3fs\n',toc);

%% Build Geometry for the Single Wall Simulation
fprintf('Building Thermal Model... ');
tic;
thermalmodel = createpde('thermal','transient');
fprintf('%1.3fs\n',toc);

fprintf('Building PDE Mesh Nodes... ');
tic;
Base = [3 4 0 base_length base_length 0 0 0 base_thick base_thick]';

%Preallocate Array Size for Layers
layers = zeros(10, num_layer.*nodeperlayer);
n = 1;

%Index Array Size Based on Node Numbers
for i = 1:num_layer
    bot_layerheight = ((i-1)*(node_thick))+base_thick;
    x_offset = (base_length - wall_length)./2;
    
    for j = 0:nodeperlayer-1
        layers(:,n) = [  3
                                4
                                x_offset+(j.*node_length)
                                x_offset+((1+j).*node_length)
                                x_offset+((1+j).*node_length)
                                x_offset+(j.*node_length)
                                bot_layerheight
                                bot_layerheight
                                bot_layerheight+node_thick
                                bot_layerheight+node_thick];
        
        n = n+1;
    end%for j

end%for i

fprintf('%1.3fs\n',toc);

fprintf('Generating PDE Mesh... ');
% Assemble
gd = [Base layers];
dl = decsg(gd);
mesh = geometryFromEdges(thermalmodel,dl);
m = generateMesh(thermalmodel,'Hmin',20);
return;
fprintf('%1.3fs\n',toc);
%% Visualize Edges

% fprintf('Plotting Faces... ');
% tic;
% figure
% pdegplot(mesh,'FaceLabels','on')
% grid on
% xlim([-10 320]);
% ylim([-10 100]);
% title 'Single Node Geometry Check';
% fprintf('%1.3fs\n',toc);

%% Build Node Properties and Assign them Activation times
% Build Array containing the time for activation for each node - this
% probably needs to be made smarter so that I don't fuck it up later

fprintf('Setting all nodes to melt temperature... ');
tic;
% Pre allocate all Nodes other than the baseplate to melting temperature
for i = 2:mesh.NumFaces
    thermalIC(thermalmodel,melt_temp,'Face',i);
end
fprintf('%1.3fs\n',toc);

fprintf('Calculating Node Activation Times... ');
tic;
% Parameter Calculation
global kActivateTime;
kActivateTime = zeros(1, num_layer.*nodeperlayer);
% Skip Face 1 because it is baseplate and initial condition of room temp...
n = 2; 
%initialize time for 1 second to give a 0 starting condition
time = 1;

for i = 0:num_layer-1
    t_offset = layerwait;
    
    for j = 1:nodeperlayer
        timepernode = (node_length./TS);
        if j == 1
            time = t_offset+time;
        else
            time = time+timepernode;
        end
        kActivateTime(n) = time;
        n = n+1;
    end
end
kActivateTime(2:end) = kActivateTime(2:end) - layerwait;

fprintf('%1.3fs\n',toc);

%%
fprintf('Setting Node Properties... ');
tic;

thermalProperties(thermalmodel,'ThermalConductivity', k,'MassDensity', density, ...
                                'SpecificHeat', Cp);
thermalIC(thermalmodel,ambient_T,'Face',1);

for i = 2:length(kActivateTime)
    thermalProperties(thermalmodel,'Face',i,'ThermalConductivity',@NodeActivationFunction, ...
        'MassDensity', density, 'SpecificHeat', Cp);
end%for i

fprintf('%1.3fs\n',toc);

%% Rerun Simulation
fprintf('Running Simulation... ');
tic;

t_sim = 60;
sim_fps = 12;
tlist = linspace(0,t_sim,t_sim*sim_fps);
% tlist = 0:0.01:30;
R = solve(thermalmodel,tlist);
T = R.Temperature;

fprintf('%1.3fs\n',toc);
return;
%% Plot Images for Visual Check

figure
hold on
pdeplot(thermalmodel,'XYData', T(:,size(tlist)),'ColorMap','hot');
xlim([-50 350]);
ylim([-50 125]);

title 'Melting Temperature Initial Condition lump end Second';
hold off

%% Multi-Image Playing

tic
disp(' ')
disp('Playing a video of the frames and window given as input.')
disp('Press "CTRL + C" to end video.')

for n=30:1:tlist(end)
    %size(frames,3)
    hold on
    pdeplot(thermalmodel,'XYData',T(:,n),'ColorMap','autumn');
    xlim([-20 350]);
    ylim([0 125]);

    title(sprintf('PDEPlot Melting Temperature Time: %s',num2str(n)))
    pause(0.01)     
end

toc


%% Check Properties and See if Assigned Correctly

mpaFace1 = findThermalProperties(thermalmodel.MaterialProperties,'Face',1)
mpaFace2 = findThermalProperties(thermalmodel.MaterialProperties,'Face',2)

%% Visualze Temp

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 80,  15);

figure(102)
hold on
plot(tlist, T(nid,:)); 

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 89,  15);
plot(tlist, T(nid,:)); 

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 96,  15);
plot(tlist, T(nid,:)); 

grid on
title('Node Locations Temperature Plot over Time for Sim Times');
xlabel('Time of Simulation (Max 40 seconds)');
ylabel('Temperature (Kelvin)');
legend('node1','node2','node3') 
xlim([-.1 tlist(end)]);
hold off

%% Visualize Temperature at Each Node for Full Time Period

getClosestNode = @(p,x,y) min((p(1,:) - x).^2 + (p(2,:) - y).^2);
[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 152.4,  13.825);

figure(102)
hold on
plot(tlist, T(nid,:)); 

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 152.4, 6);
plot(tlist, T(nid,:)); 

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 32, 6);
plot(tlist, T(nid,:)); 

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 270, 6);
plot(tlist, T(nid,:));

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 152.4, 50);
plot(tlist, T(nid,:));

[~,nid] = getClosestNode( thermalmodel.Mesh.Nodes, 152.4, 78);
plot(tlist, T(nid,:));

grid on
title('Node Locations Temperature Plot over Time for Full Time');
xlabel('Time of Simulation (Max 1 Second)');
ylabel('Temperature (Kelvin)');
legend('Center of Node','Center of Plate','Left Center Plate','Right Center Plate','Middle of Wall','Top of Wall') 
xlim([-.1 1]);
hold off
