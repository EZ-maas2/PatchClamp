% Elizaveta Zavialova for Neurobiology of isolated tissue practical 2026
%% Load in the timeseries
[p,d] = uigetfile();
ts = load(fullfile(d,p)); % Y entries contain the data abnout the current and volatge change
%% Constants
fs = string(fieldnames(ts));
fs = fs(contains(fs, "Y"));
NUM_SIGNALS = numel(fs);

sf = 20000; %20kHz 
currentStart = 2000;
currentStop = 4000;



%% Retrieve currents and volatges
t_ms = ts.T1 * 1000;

current_traces = nan(NUM_SIGNALS, numel(t_ms));
og_currents = nan(NUM_SIGNALS, 1);
voltages = nan(NUM_SIGNALS, numel(t_ms));

for f = 1:NUM_SIGNALS
    rec = ts.(fs(f));
    current_traces(f, :) = rec(:, 1);
    og_currents(f) = mean(current_traces(f, currentStart:currentStop)); 
    voltages(f, :) = rec(:, 2);
end
%% Current correction term
% Find the internal current 
if false
    stds = std(voltages, 0, 2);
    neutral_current_ix = find(stds < 2);
    compensatory_current = og_currents(neutral_current_ix);
    
    currents = og_currents - compensatory_current;
    current_traces = current_traces - compensatory_current;
end

%%
% Show the profile of the stimulation
volt_f = figure;
curr_f = figure;

cmap = colormap(parula(15));

for f = 1:NUM_SIGNALS
    
 
    %voltage
    figure(volt_f)
    plot(t_ms, voltages(f, :), 'Color', cmap(f, :), "DisplayName", sprintf("%1.2f nA",currents(f))); hold on; 

    %current
    figure(curr_f)
    plot(t_ms, current_traces(f, :), 'Color', cmap(f, :)); hold on; 
    plot([200, 200], [- 0.3, 0.3], 'k--'); hold on;
    plot([100, 100], [- 0.3, 0.3], 'k--'); hold on;



end

ylabel("Current, nA");
xlabel("Time, ms"); 

figure(volt_f)
ylabel("Voltage, mV");
xlabel("Time, ms"); 
legend('Location','northeastoutside');

%% Save figures
figure(curr_f);
exportgraphics(gcf, "Current.png");
figure(volt_f);
exportgraphics(gcf, replace(p, ".mat", "")+"_Voltage.png");


%% Find the baseline values
baselines = mean(voltages(:, 1:currentStart), 2); % mean of rows, cause every row is a timeseries
smoothed_voltages = movmean(voltages, 100, 2); % 400 datapoints window = 20 ms smoothing


%% Calculate the AP amplitude & threshold

% Find the traces with APs
ap_mask = false(NUM_SIGNALS, 1);
for i = 1:NUM_SIGNALS
    ap_mask(i) = sum(smoothed_voltages(i, :) >=  baselines(i) + 30) > 0;
end
ap_ixs = find(ap_mask);

% Out of all of them, select the one with the minimal current
rheobase = min(currents(ap_mask)); % min current necessary to elicit a spike
rheobase_ix = find(currents == rheobase);

% Threshold - sharp elbow where we go from slwoly charging to spiking
[~,biggest_der_ix] = max(diff(diff(smoothed_voltages(rheobase_ix, :)))); % biggest second derivative
%% Count the number of spikes
zscored = zscore(voltages(:, currentStart:currentStop), 0,2); % look at how zscores represent the aps
zscored_f = figure;
num_spikes = nan(NUM_SIGNALS, 1);
isi = nan(NUM_SIGNALS, 1);

for i = 1:numel(ap_ixs)
    figure;
    f = ap_ixs(i);
    title(sprintf("Stimulation current is %1.2f nA",currents(f)));

    lmax = islocalmax(zscored(f, :)) & zscored(f,:) > 1;
    % Calculate interspike interval (ms)
    spike_ixs = find(lmax);
    isi(f) = 1000*mean(diff(spike_ixs))/sf;
    num_spikes(f) = sum(lmax);
    plot(t_ms(currentStart:currentStop), zscored(f, :), 'Color', cmap(f, :), "DisplayName", sprintf("%1.2f nA",currents(f))); hold on; 
    
    maxes = zscored(f, :);
    maxes(~lmax) = nan;
    plot(t_ms(currentStart:currentStop),maxes, "r*"); hold off;
    xlabel("Time, ms");
end


%% Plot the threshold on the voltage figure
figure(volt_f);
Threshold_voltage = voltages(rheobase_ix, biggest_der_ix);
plot(biggest_der_ix / sf * 1000, Threshold_voltage, "DisplayName", "Threshold", "Marker","x", 'MarkerSize',8, 'Color','Black'); hold on;

%% Calculate the firing rate
% Interspike interval


%% Calculate the tau (time constant) via finding 2/3 of the exponential decay
% tau = time at which two thirds of the signal has been lost

[~, max_ix] = max(abs(smoothed_voltages - baselines), [], 2); % find teh biggest deviation from the baseline
max_v = nan(NUM_SIGNALS, 1);
two_thirds = nan(NUM_SIGNALS, 1);
tau = nan(NUM_SIGNALS, 1);
intersection_ixs = nan(NUM_SIGNALS, 1);

for ix = 1:NUM_SIGNALS
    max_v(ix) = smoothed_voltages(ix, max_ix(ix));
    two_thirds(ix) = (max_v(ix) - baselines(ix)) * 2/3 + baselines(ix);

    % find the intersection of voltages and the two_thirds
    [~, intersection_ix] = min(abs(smoothed_voltages(ix, :) - two_thirds(ix)));
    intersection_ixs(ix) = intersection_ix;
    tau(ix) = intersection_ix / sf; % convert tau to seconds

end



%% Calculate the Membrane resistance R based on the segment of the voltage trace where dV/dt = 0 (all current is through R)

stimulationSegment = smoothed_voltages(:, currentStart:currentStop);

% take a row-wise temporal derivative - find dV/dt approx 0
mask = [abs(diff(stimulationSegment, 1, 2)) < 0.01 false(NUM_SIGNALS, 1)]; % padded with logical 0 

% only take last 25 ms of the stimulation window
mask(1:round(3/4*numel(mask))) = false;

% calculate R as mean voltage in places where dV/dt = 0 divided by
% stimulation current
R = mean(stimulationSegment .* mask, 2) ./ currents;


% Plot the voltage traces with the regions that 
if true
    mask_f = figure();
    for i = 1:NUM_SIGNALS
        masked = stimulationSegment(i, :);
        masked(mask(i,:)) = nan;
        plot(masked); hold on;
    end
end
%% Calculate C as tau/R
C = tau./R;
%% Plot voltage traces, regions taken as dV/dt, intersections for the tau calculation and AP