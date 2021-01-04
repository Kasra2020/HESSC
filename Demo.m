clear all;
close all;
clc;

%% Add the dependencies
addpath (genpath('./Functions'))

%% Automatic HESSC
depth = 3; % The depth for exploring the data for the automatic version of HESSC.
[CSmap,Ptime, n_cluster] = HESSC_fauto(img,depth);

%% Manual HESSC
% n_cluster = 5; % Numeber of clusters for the manual version of HESSC.
% [CSmap, Ptime] = HESSC_fm(img,n_cluster);
