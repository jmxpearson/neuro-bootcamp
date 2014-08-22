%% load it up
cd('~/code/neuro-bootcamp-2014')
addpath('~/code/electrophysiology')
load sample_spikes
load sample_lfp
load sample_events
whos
%% check some data
sr
numel(lfp) / sr
numel(lfp) / (sr * 60)

mean(lfp)
median(lfp)

%% take a look
plot(lfp)
plot(lfp(1:100))

plot(lfp, 'k-', 'linewidth', 1.5)

lfplen = length(lfp);
taxis = (1:lfplen) / sr;  % what happens if we omit parens?
plot(taxis, lfp, 'k')
xlabel('Time (s)')
ylabel('Voltage (mV)')

xlim([0 20])

plot([0 20], [0.2 0.2], 'r') %whoops!

%% try again
figure
hold all
lfplen = length(lfp);
taxis = (1:lfplen) / sr;  % what happens if we omit parens?
plot(taxis, lfp, 'k')
xlabel('Time (s)')
ylabel('Voltage (mV)')

xlim([0 20])

plot([0 20], [0.2 0.2], 'r')
hold off

%% how would we make a psth from spike times?

% one approach: bin spikes, then cut
% alternative: make relative time bins, then histogram

startT = -1.5;
stopT = 1.0;
binsize = 0.001;  % 10 ms bins

taxis = startT:binsize:stopT;

for ind = 1:numel(events)
   this_taxis = events(ind) + taxis;  % moveable window
   raster(ind, :) = histc(times, this_taxis);
end

%% plot
colormap('gray')
imagesc(raster)

%% how do we make a psth?
% histogram: sum across trials; trials are rows
psth = sum(raster); % same as sum(raster, 1)

bar(taxis, psth, 'k')

%% or try a plot
plot(taxis, psth)

smooth_time = 0.050;
smooth_bins = smooth_time / binsize;
plot(taxis, smooth(psth, smooth_bins))  % note the edge artifacts!
ylim([10, 20])