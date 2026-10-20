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
currents = nan(NUM_SIGNALS, 1);
voltages = nan(NUM_SIGNALS, numel(t_ms));

for f = 1:NUM_SIGNALS
    rec = ts.(fs(f));
    current_traces(f, :) = rec(:, 1);
    currents(f) = mean(current_traces(f, currentStart:currentStop)); 
    voltages(f, :) = rec(:, 2);
end
%% Current correction term
% If the baseline was not zero, correct it to zero
currents = currents - mean(current_traces(1:currentStart - 10));

%% Color map
cmap = colormap(parula(16));
%%
% Show the profile of the stimulation
volt_f = figure;
curr_f = figure;


ls  = {'--','-.','-','--','-.','-', '--','-.','-', '--','-.','-', '--','-.','-'};
mk  = {'none','none','none','none','none','none', 'none','none','none','none','none','none','none','none','none'};%,'o','s','^', 'o','s','^', 'x','+','d'};


for f = 1:NUM_SIGNALS
    
    %voltage
    figure(volt_f)
    plot(t_ms, voltages(f, :), 'LineStyle',ls{f}, 'Marker',mk{f}, 'Color', cmap(f, :), "DisplayName", sprintf("ΔI = %1.2f nA",currents(f))); hold on; 

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
exportgraphics(gcf, replace(p, ".mat", "")+"_Current.png");
figure(volt_f);
exportgraphics(gcf, replace(p, ".mat", "")+"_Voltage.png");


%% Find the baseline values
baselines = mean(voltages(:, 1:currentStart), 2); % mean of rows, cause every row is a timeseries
smoothed_voltages = movmean(voltages, 10, 2); % 400 datapoints window = 0.5 ms smoothing
stimulationSegment = smoothed_voltages(:, currentStart:currentStop);

Baseline.avg = mean(baselines);
Baseline.std = std(baselines);

%% Calculate the AP amplitude & threshold
% Find the traces with APs
zscored = zscore(voltages(:, currentStart:currentStop), 0, 2); % look at how zscores represent the aps
ap_mask = false(NUM_SIGNALS, 1);
neutral_mask = false(NUM_SIGNALS, 1);
for i = 1:NUM_SIGNALS
    if std(voltages(i, :)) > 2 % if this is a non-neutral trace
        ap_mask(i) = sum(zscored(i, :) > 1.5 & voltages(i, currentStart:currentStop) > baselines(i)) > 0; % if there are spikes
    else
        neutral_mask(i) = true;
    end
end

ap_ixs = find(ap_mask);

% Out of all of them, select the one with the minimal current
rheobase = min(currents(ap_mask)); % min current necessary to elicit a spike
rheobase_ix = find(currents == rheobase);
%% Count the number of spikes

num_spikes = nan(NUM_SIGNALS, 1);
isi = nan(NUM_SIGNALS, 1);

zscored_f = figure('Position', [100 100 1400 900],'WindowState', 'maximized');
t = tiledlayout(numel(ap_ixs), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
spike_ixs = cell(NUM_SIGNALS, 1); 

for i = 1:numel(ap_ixs)
    nexttile;
    f = ap_ixs(i);

    lmax = islocalmax(zscored(f, :)) & zscored(f,:) > 1;

    % Calculate interspike interval (ms)
    spike_ixs{f} = find(lmax);
    isi(f) = 1000*mean(diff(spike_ixs{f}))/sf;


    num_spikes(f) = sum(lmax);
    plot(t_ms(currentStart:currentStop), voltages(f, currentStart:currentStop), 'Color', cmap(f, :), "DisplayName",  sprintf("ΔI = %1.2f nA",currents(f))); hold on; 
    
    maxes = voltages(f, currentStart:currentStop);
    maxes(~lmax) = nan;
    plot(t_ms(currentStart:currentStop),maxes, "k*"); hold off;
    xlabel("Time, ms");
    ylabel("Voltage, mV")
    title(sprintf("ΔI = %1.2f nA",currents(f)));
    spike_ixs{f} = spike_ixs{f} + currentStart;

end
exportgraphics(gcf, replace(p, ".mat", "")+"_Spikes.png")

%% Threshold - sharp elbow where we go from slowly charging to spiking

% Do it for each spike
spikes_in_rheobase = spike_ixs{rheobase_ix};
threshold = nan(numel(spikes_in_rheobase), 1);
biggest_der_ix = nan(numel(spikes_in_rheobase), 1);


for ix = 1:numel(spikes_in_rheobase)
    spike_ix = spikes_in_rheobase(ix); 
    % search from 10ms before AP up to AP
    ms2 = 2 * sf /1000;
    window = voltages(rheobase_ix, spike_ix-ms2:spike_ix);
    [~, bdix] = max(diff(diff(window)));
    bdix = bdix - 5; % go 0.25 ms before, because second derivative underestimates 
    threshold(ix) = window(bdix);
    biggest_der_ix(ix) = spike_ix-ms2 + bdix - 1;
    figure;
    plot(window);hold on;
    plot(bdix, window(bdix), "r*");
end

Threshold = min(threshold);
%% Plot the threshold on the voltage figure
figure(volt_f);
for i = 1:numel(biggest_der_ix)
    tv = voltages(rheobase_ix, biggest_der_ix(i));
    plot(t_ms(biggest_der_ix(i)), tv, "DisplayName", "Threshold", "Marker","x", 'MarkerSize',12, 'Color', 'Black');hold on;
end

%% Calculate the tau (time constant) via finding 1-1/e of the exponential decay
% tau = time at which two thirds of the signal has been lost

[~, max_ix] = max(abs(stimulationSegment - baselines), [], 2); % find teh biggest deviation from the baseline

max_v = nan(NUM_SIGNALS, 1);
two_thirds = nan(NUM_SIGNALS, 1);
tau = nan(NUM_SIGNALS, 1);
intersection_ixs = nan(NUM_SIGNALS, 1);

for ix = 1:NUM_SIGNALS
    max_v(ix) = stimulationSegment(ix, max_ix(ix));
    two_thirds(ix) = (max_v(ix) - baselines(ix)) * (1-exp(-1)) + baselines(ix);

    % find the intersection of voltages and the two_thirds
    [~, intersection_ix] = min(abs(stimulationSegment(ix, :) - two_thirds(ix)));
    intersection_ixs(ix) = intersection_ix;

    figure(volt_f);
    plot(t_ms(intersection_ix+currentStart), stimulationSegment(ix, intersection_ix), "k*", "DisplayName","Intersection"); hold on;

    tau(ix) = intersection_ix/sf; % convert tau to seconds

end



%% Calculate the Membrane resistance R based on the segment of the voltage trace where dV/dt = 0 (all current is through R)


% take a row-wise temporal derivative of hyperpolarized segments - find dV/dt approx 0
relevant_traces = ~neutral_mask & ~ap_mask;

mask = [abs(diff(stimulationSegment, 1, 2)) < 0.01 false(NUM_SIGNALS, 1)]; % padded with logical 0 


% only take last 25 ms of the stimulation window
mask(1:round(3/4*numel(mask))) = false;
mask(~relevant_traces) = false;

dV = nan(NUM_SIGNALS, 1);


% Plot the voltage traces with the regions that have dV/dt = 0
time = t_ms(currentStart:currentStop);
if true
    mask_f = figure;
    for i = 1:NUM_SIGNALS
        masked = stimulationSegment(i, :);
        masked(mask(i,:)) = nan;

        plot(time, masked, 'Color',cmap(i, :)); hold on;
        plot(time(mask(i,:)), stimulationSegment(i, mask(i,:)), "r-.");
        
        % take voltages from dV/dt = 0 region
        dV(i) = mean(stimulationSegment(i, mask(i,:)), 2);
    end
end

% calculate R as mean voltage in places where dV/dt = 0 divided by
% stimulation current
dV = dV - Baseline.avg;
R = dV./ currents; % units are MOhm (10^6)

%% Plot steady-state current-voltage relationship
vi_fig = figure;
title("I-V relationship");
plot(currents, dV);
xlabel("ΔI, nA");
ylabel("ΔV, mV");
exportgraphics(gcf, replace(p, ".mat", "")+"_IV.png")


%% Calculate C as tau/R
C = tau./(R * 1e6); % R is in MOhm
%% Combined Current-Voltage plot

final_figure = figure('Position', [100 100 1400 900],'WindowState', 'maximized');
title("Overview of the responses to the step currents");

for f = 1:NUM_SIGNALS
    %f = NUM_SIGNALS - s + 1;
    %voltage
    subplot(1, 5, [1, 3]);
    plot(t_ms, voltages(f, :), 'LineStyle',ls{f}, 'Marker',mk{f}, 'Color', cmap(f, :), "DisplayName", sprintf("ΔI = %1.2f nA",currents(f))); hold on; 
   

    %current
    subplot(1, 5, [4,5]);
    plot(t_ms, current_traces(f, :), 'Color', cmap(f, :)); hold on; 
    plot([200, 200], [- 0.3, 0.3], 'k--'); hold on;
    plot([100, 100], [- 0.3, 0.3], 'k--'); hold on;


end

subplot(1, 5, [1, 3]);
title("Voltage traces",'FontSize', 23);
l = legend('Location','northeastoutside', 'FontSize', 12);
ylabel("Voltage, mV",'FontSize', 18);
xlabel("Time, ms",'FontSize', 18);

subplot(1, 5, [4,5]);
title("Step currents",'FontSize', 23);
ylabel("Current, nA",'FontSize', 18);
xlabel("Time, ms",'FontSize', 18); 




