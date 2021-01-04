clear all;
close all;
clc;

%% Add the dependencies
addpath (genpath('./Functions'))

%% Setting the essential parameters for each function
depth = 3; % The depth for exploring the data for the automatic version of HESSC.
n_cluster = 5; % Numeber of clusters for the manual version of HESSC.

%% Automatic HESSC
% CSmap = HESSC_fauto(img,depth);

%% Manual 
CSmap = HESSC_fm(img,5);

figure;imagesc(CSmap);axis off;axis image;title('HESSC')