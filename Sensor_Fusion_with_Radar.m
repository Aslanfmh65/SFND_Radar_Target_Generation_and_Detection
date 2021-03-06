clear all 
clc;

%% Define an empty scenario
% Create scenario object
scenario = drivingScenario;
% scenario will generate/update at every 0.01s
scenario.SampleTime = 0.01;

%% Create road
% Create a highway of 500 meters length, 7.5 meter width, with two lane
roadCenters = [0 0; 50 0; 100 0; 250 20; 500 40];
roadWidth = 7.2; % Two lanes, each 3.6 meters
road(scenario, roadCenters, roadWidth);

%% Create vehicle 
% Create the ego vehicle that travels at 25 m/s along the road.
egoCar = vehicle(scenario, 'ClassID', 1);
path(egoCar, roadCenters(2:end,:) - [0 1.8], 25); % On right lane

% Add a car in front of the ego vehicle.
leadCar = vehicle(scenario, 'ClassID', 1);
path(leadCar, [70 0; roadCenters(3:end,:)] - [0 1.8], 25); % On right lane

% Add a car that travels at 35 m/s along the road and passes ego vehicle.
passingCar = vehicle(scenario, 'ClassID', 1);
waypoints = [0 -1.8; 50 1.8; 100 1.8; 250 21.8; 400 32.2; 500 38.2];
path(passingCar, waypoints, 35);

% Add a car behind the ego vehicle
chaseCar = vehicle(scenario, 'ClassID', 1);
path(chaseCar, [25 0; roadCenters(2:end,:)] - [0 1.8], 25); % On right lane

%% Define Radar
sensors = cell(2,1);

% Front-facing long-range radar sensor at car's front bumper
sensors{1} = radarDetectionGenerator(...
          'SensorIndex', 1, ...
          'Height', 0.2, ...
          'MaxRange', 174, ...
          'SensorLocation',[egoCar.Wheelbase + egoCar.FrontOverhang, 0],...
          'FieldOfView', [20, 5]);
  
% Rear-facing long-range radar sensor at car's rear bumper
sensors{2} = radarDetectionGenerator(...
          'SensorIndex', 2, ...
          'Height', 0.2, ...
          'Yaw', 180, ...
          'SensorLocation',[-egoCar.RearOverhang, 0],...
          'MaxRange', 174,...
          'FieldOfView', [20, 5]);
          
% Rear-left-facing short-range radar sensor at the left rear wheel of car
sensors{3} = radarDetectionGenerator(...
          'SensorIndex', 3, ...
          'Height', 0.2, ...
          'Yaw', 120, ...
          'SensorLocation',[0, egoCar.Width/2],...
          'MaxRange', 30,...
          'ReferenceRange', 50,...
          'FieldOfView', [90, 5],...
          'AzimuthResolution', 10,...
          'RangeResolution', 1.25);

% Rear-right-facing short-range radar sensor at the right rear wheel of car
sensors{4} = radarDetectionGenerator(...
          'SensorIndex', 4, ...
          'Height', 0.2, ...
          'Yaw', -120, ...
          'SensorLocation',[0, -egoCar.Width/2],...
          'MaxRange', 30,...
          'ReferenceRange', 50,...
          'FieldOfView', [90, 5],...
          'AzimuthResolution', 10,...
          'RangeResolution', 1.25);

% Front-left-facing short-range radar sensor at the left front wheel of car
sensors{5} = radarDetectionGenerator(...
          'SensorIndex', 5, ...
          'Height', 0.2, ...
          'Yaw', 60, ...
          'SensorLocation',[egoCar.Wheelbase, egoCar.Width/2],...
          'MaxRange', 30,...
          'ReferenceRange', 50,...
          'FieldOfView', [90, 5],...
          'AzimuthResolution', 10,...
          'RangeResolution', 1.25);

% Front-right-facing short-range radar sensor at the right front wheel of
% car
sensors{6} = radarDetectionGenerator(...
          'SensorIndex', 6, ...
          'Height', 0.2, ...
          'Yaw', -60, ...
          'SensorLocation',[egoCar.Wheelbase, -egoCar.Width/2],...
          'MaxRange', 30,...
          'ReferenceRange', 50,...
          'FieldOfView', [90, 5],...
          'AzimuthResolution', 10,...
          'RangeResolution', 1.25);
      
% Front-facing camera located at front windshield.
sensors{7} = visionDetectionGenerator(...
          'SensorIndex', 7,...
          'FalsePositivesPerImage', 0.1, ...
          'SensorLocation', [0.75*egoCar.Wheelbase 0],...
          'Height', 1.1);

% Rear-facing camera located at rear windshield.
sensors{8} = visionDetectionGenerator(...
          'SensorIndex', 8,...
          'FalsePositivesPerImage', 0.1, ...
          'SensorLocation', [0.2*egoCar.Wheelbase 0],...
          'Height', 1.1, 'Yaw', 180);

profiles = actorProfiles(scenario);
for m = 1:numel(sensors)
    sensors{m}.ActorProfiles = profiles;
end
         
%% Create MultiObejectTracker
% Create a multiObjectTracker to track the vehicles that are close to 
% the ego vehicle
tracker = multiObjectTracker(...
          'FilterInitializationFcn', @initSimDemoFilter, ...
          'AssignmentThreshold', 30,...
          'ConfirmationParameters', [4 5],...
          'NumCoastingUpdates', 5);

positionSelector = [1 0 0 0; 0 0 1 0]; % Position selector
velocitySelector = [0 1 0 0; 0 0 0 1]; % Velocity selector

% Create the display and return a handle to the bird's-eye plot
BEP = createDemoDisplay(egoCar, sensors);

%% Simulate Scenario
toSnap = true;
while advance(scenario) && ishghandle(BEP.Parent)    
    % Get the scenario time
    time = scenario.SimulationTime;

    % Get the position of the other vehicle in ego vehicle coordinates
    ta = targetPoses(egoCar);

    % Simulate the sensors
    detections = {};
    isValidTime = false(1,6);
    for i = 1:6
        [sensorDets,numValidDets,isValidTime(i)] = sensors{i}(ta, time);
        if numValidDets
            for j = 1:numValidDets
                % Vision detections do not report SNR. The tracker requires
                % that they have the same object attributes as the radar
                % detections. This adds the SNR object attribute to vision
                % detections and sets it to a NaN.
                if ~isfield(sensorDets{j}.ObjectAttributes{1}, 'SNR')
                    sensorDets{j}.ObjectAttributes{1}.SNR = NaN;
                end
            end
            detections = [detections; sensorDets]; %#ok<AGROW>
        end
    end

    % Update the tracker if there are new detections
    if any(isValidTime)
        vehicleLength = sensors{1}.ActorProfiles.Length;
        detectionClusters = clusterDetections(detections, vehicleLength);
        confirmedTracks = updateTracks(tracker, detectionClusters, time);

        % Update bird's-eye plot
        updateBEP(BEP, egoCar, detections, confirmedTracks, ...
                  positionSelector, velocitySelector);
    end

    % Snap a figure for the document when the car passes the ego vehicle
    if ta(1).Position(1) > 0 && toSnap
        toSnap = false;
        snapnow
    end
end

%% Kalman Filter Function
function filter = initSimDemoFilter(detection)

% Use a 2-D constant velocity model to initialize a trackingKF filter.
% The state vector is [x;vx;y;vy]
% The detection measurement vector is [x;y;vx;vy]
% As a result, the measurement model is H = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1]

H = [1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1];
filter = trackingKF('MotionModel', '2D Constant Velocity', ...
 'State', H' * detection.Measurement, ...
 'MeasurementModel', H, ...
 'StateCovariance', H' * detection.MeasurementNoise * H, ...
 'MeasurementNoise', detection.MeasurementNoise);
end

%% Clustering Detection Function
function detectionClusters = clusterDetections(detections, vehicleSize)
N = numel(detections);
distances = zeros(N);

for i = 1:N
    for j = i+1:N
        if detections{i}.SensorIndex == detections{j}.SensorIndex
            distances(i,j) = norm(detections{i}.Measurement(1:2) - ...
                                  detections{j}.Measurement(1:2));
        else
            distances(i,j) = inf;
        end
    end
end

leftToCheck = 1:N;
i=0;
detectionClusters = cell(N,1);
while ~isempty(leftToCheck)
    underConsideration = leftToCheck(1);
    clusterInds = (distances(underConsideration, leftToCheck) < vehicleSize);
    detInds = leftToCheck(clusterInds);
    clusterDets = [detections{detInds}];
    clusterMeas = [clusterDets.Measurement];
    meas = mean(clusterMeas,2);
    meas2D = [meas(1:2); meas(4:5)];
    i = i+1;
    detectionClusters{i} = detections{detInds(1)};
    detectionClusters{i}.Measurement  = meas2D;
    leftToCheck(clusterInds) = [];
    
end
   
detectionClusters(i+1:end) = [];

% Since the detections are now for clusters, modify the noise to represent
% that they are of the whole car
for i = 1:numel(detectionClusters)
    measNoise(1:2,1:2) = vehicleSize^2 * eye(2);
    measNoise(3:4,3:4) = eye(2) * 100 * vehicleSize^2;
    detectionClusters{i}.MeasurementNoise = measNoise;
end
end

%% Display Function
function BEP = createDemoDisplay(egoCar, sensors)
    % Make a figure
    hFigure = figure('Position', [0, 0, 1200, 640], 'Name',...
                     'Sensor Fusion with Synthetic Data Example');
    % Moves the figure to the left and a little down from the top
    movegui(hFigure, [0 -1]); 

    % Add a car plot that follows the ego vehicle from behind
    hCarViewPanel = uipanel(hFigure, 'Position', [0 0 0.5 0.5],...
                           'Title', 'Chase Camera View');
    hCarPlot = axes(hCarViewPanel);
    chasePlot(egoCar, 'Centerline', 'on', 'Parent', hCarPlot);

    % Add a car plot that follows the ego vehicle from a top view
    hTopViewPanel = uipanel(hFigure, 'Position', [0 0.5 0.5 0.5],...
                            'Title', 'Top View');
    hCarPlot = axes(hTopViewPanel);
    chasePlot(egoCar, 'Centerline', 'on', 'Parent', hCarPlot,...
              'ViewHeight', 130, 'ViewLocation', [0 0], 'ViewPitch', 90);

    % Add a panel for a bird's-eye plot
    hBEVPanel = uipanel(hFigure, 'Position', [0.5 0 0.5 1],...
                       'Title', 'Bird''s-Eye Plot');

    % Create bird's-eye plot for the ego car and sensor coverage
    hBEVPlot = axes(hBEVPanel);
    frontBackLim = 60;
    BEP = birdsEyePlot('Parent', hBEVPlot, 'Xlimits',...
                      [-frontBackLim frontBackLim], 'Ylimits', [-35 35]);

    % Plot the coverage areas for radars
    for i = 1:6
        cap = coverageAreaPlotter(BEP,'FaceColor','red','EdgeColor','red');
        plotCoverageArea(cap, sensors{i}.SensorLocation,...
            sensors{i}.MaxRange, sensors{i}.Yaw, sensors{i}.FieldOfView(1));
    end
    
    % Plot the coverage areas for vision sensors
    for i = 7:8
        cap = coverageAreaPlotter(BEP,'FaceColor','blue','EdgeColor','blue');
        plotCoverageArea(cap, sensors{i}.SensorLocation,...
            sensors{i}.MaxRange, sensors{i}.Yaw, 45);
    end

    % Create a vision detection plotter put it in a struct for future use
    detectionPlotter(BEP, 'DisplayName','vision', 'MarkerEdgeColor',...
                    'blue', 'Marker','^');

    % Combine all radar detctions into one entry and store it for later update
    detectionPlotter(BEP, 'DisplayName','radar', 'MarkerEdgeColor','red');

    % Add road borders to plot
    laneBoundaryPlotter(BEP, 'DisplayName','road', 'Color', [.75 .75 0]);

    % Add the tracks to the bird's-eye plot. Show last 10 track updates.
    trackPlotter(BEP, 'DisplayName','track', 'HistoryDepth',10);

    axis(BEP.Parent, 'equal');
    xlim(BEP.Parent, [-frontBackLim frontBackLim]);
    ylim(BEP.Parent, [-40 40]);

    % Add an outline plotter for ground truth
    outlinePlotter(BEP, 'Tag', 'Ground truth');
end



function updateBEP(BEP, egoCar, detections, confirmedTracks, psel, vsel)
    % Update road boundaries and their display
    rb = roadBoundaries(egoCar);
    plotLaneBoundary(findPlotter(BEP,'DisplayName','road'),rb);

    % update ground truth data
    [position, yaw, length, width, originOffset, color] = targetOutlines(egoCar);
    plotOutline(findPlotter(BEP,'Tag','Ground truth'), position, yaw,...
             length, width, 'OriginOffset', originOffset, 'Color', color);

    % Prepare and update detections display
    N = numel(detections);
    detPos = zeros(N,2);
    isRadar = true(N,1);
    for i = 1:N
        detPos(i,:) = detections{i}.Measurement(1:2)';
        if detections{i}.SensorIndex > 6 % Vision detections
            isRadar(i) = false;
        end
    end
    plotDetection(findPlotter(BEP,'DisplayName','vision'), detPos(~isRadar,:));
    plotDetection(findPlotter(BEP,'DisplayName','radar'), detPos(isRadar,:));

    % Prepare and update tracks display
    trackIDs = {confirmedTracks.TrackID};
    labels = cellfun(@num2str, trackIDs, 'UniformOutput', false);
    [tracksPos, tracksCov] = getTrackPositions(confirmedTracks, psel);
    tracksVel = getTrackVelocities(confirmedTracks, vsel);
    plotTrack(findPlotter(BEP,'DisplayName','track'), tracksPos, ...
                          tracksVel, tracksCov, labels);
end