% Elizaveta Zavialova for Neurobiology of isolated tissue practical 2026
%% Load in the timeseries
[p,d] = uigetfile();
ts = load(fullfile(d,p)); % Y entries contain the data abnout the current and volatge change
%%
sf = 20000; %20kHz 
currentStart = 2000;
%%
% Show the profile of the stimulation
volt_f = figure;
curr_f = figure;
fs = string(fieldnames(ts));
fs = fs(contains(fs, "Y"));
cmap = colormap(parula(15));
t_ms = ts.T1 * 1000;
currents = nan(numel(fs), 1);
voltages = nan(numel(fs), numel(t_ms));

for f = 1:numel(fs)
    rec = ts.(fs(f));

    currents(f) = mean(rec(2000:4000, 1)); %
    voltages(f, :) = rec(:, 2);

    %voltage
    figure(volt_f)
    plot(t_ms, rec(:, 2), 'Color', cmap(f, :), "DisplayName", sprintf("%1.2f pA",currents(f))); hold on; 

    %current
    figure(curr_f)
    plot(t_ms, rec(:, 1), 'Color', cmap(f, :)); hold on; 
    plot([200, 200], [- 0.3, 0.3], 'k--'); hold on;
    plot([100, 100], [- 0.3, 0.3], 'k--'); hold on;



end

ylabel("Current, pA");
xlabel("Time, ms"); 
exportgraphics(gcf, "Current.png")

figure(volt_f)
ylabel("Voltage, mV");
xlabel("Time, ms"); 
legend('Location','northeastoutside');
exportgraphics(gcf, replace(p, ".mat", "")+"_Voltage.png")


%% Find the baseline values
baselines = mean(voltages(:, 1:currentStart), 2); % mean of rows, cause every row is a timeseries
voltages = movmean(voltages, 200, 2); % 200 datapoints window = 10 ms smoothing

% Find which voltage traces contain APs and how many APs they contain
% What do we base our jusgmnet on? Derivative (extremely big shoot up)

%% Calculate the AP amplitude


%% Calculate the firing rate
%% Calculate the tau (time constant) via fitting an exponential
% tau = time at which two thirds of the signal has lost


%% Calculate the Membrane resistance based on V-I graph

