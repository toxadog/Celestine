%% Extract starting point of each trajectory and save
clear all
cd 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.13.30 2%  spont'
load('Waves.mat')
WaveNames=fieldnames(Waves);
WavestartPoint = zeros(2,numel(WaveNames));
WavestartTime = zeros(1,numel(WaveNames));
WaveEndTime = zeros(1,numel(WaveNames));
for i=1:numel(WaveNames)
    WavestartPoint(:,i) = [Waves.(WaveNames{i}).x(1);Waves.(WaveNames{i}).y(1)];
     WavestartTime(i) =  Waves.(WaveNames{i}).t0;
     WaveEndTime(i) =  Waves.(WaveNames{i}).tn;
end

svpath = 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.13.30 2%  spont';
filename = sprintf('WavestartPoint');
strfile = fullfile(svpath, filename);
save(strfile,'WavestartPoint');

svpath = 'D:\SP_step1\Mouse 9\Experiment 2022-03-30T 18.13.30 2%  spont';
filename = sprintf('WavestartTime');
strfile = fullfile(svpath, filename);
save(strfile,'WavestartTime');

filename = sprintf('WaveEndTime');
strfile = fullfile(svpath, filename);
save(strfile,'WaveEndTime');




