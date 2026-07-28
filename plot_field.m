clear
close all
baseDir = '/scratch.global/kuma0458/c14ak1_re180/run';
c=14;
tstart=3140000000;
step = 400000;
tend=3148000000;

load(fullfile( baseDir,'grid.mat'))
%load(fullfile( baseDir,'partial_sum.mat'))
load(fullfile(baseDir,'mean_phase_fields.mat'))