%% Example : Cluster Analysis
%
% This example loads the derived model correlation matrices that were used
% in the previous SPD and Correlation paper. There are 3 matrices under the
% file 'cluster_three_corrs.mat' that are mutually distinct. In this
% example, 10 matrices are generated from each model correlation matrix by
% perturbing a little bit, leading to a total of 30 correlation matrices
% consisting of three classes. 



%% initialization
clear; close all; clc;

%  add the current directory to working path
addpath(genpath(pwd));

%  load three model correlation matrices
load("cluster_three_corrs.mat");

%% data generation
%  parameters
ncopy = 10; % number of samples per class
sd = 0.1;   % degree of perturbation

%  an empty array for recording correlation matrices
data_3d = zeros(5,5,3*ncopy);

%  create perturbed versions
for i=1:ncopy
    data_3d(:,:,i) = corraux_perturb(C1, sd);
    data_3d(:,:,i+ncopy) = corraux_perturb(C2, sd);
    data_3d(:,:,i+(2*ncopy)) = corraux_perturb(C3, sd);
end

%  "wrap" the data using 'corr_ini
data_obj = corr_init(data_3d);

%  create a true label of the data
data_lab = repelem(1:3, ncopy);

%% apply spectral clustering algorithm
%  spectral clustering have controllable parameters:
%       - K    : predefined number of clusters, and 
%       - nnbd : size of the neighborhood for bandwidth control
%  In this example, set 'nnbd=5' (which is default) and change K.

label2 = corr_specc(data_obj, "K", 2);
label3 = corr_specc(data_obj, "K", 3);
label4 = corr_specc(data_obj, "K", 4);

%% compute 2-dimensional embedding for visualization
%  we use classical MDS for its simplicity.
embed2 = corr_cmds(data_obj);

%% visualize

tiledlayout(2,2, "TileSpacing", "compact", "Padding", "tight");
nexttile
scatter(embed2(:,1), embed2(:,2), 20, data_lab, 'filled');
xlabel("DIM1"); ylabel("DIM2"); axis square; title("true label");
nexttile
scatter(embed2(:,1), embed2(:,2), 20, label2, 'filled');
xlabel("DIM1"); ylabel("DIM2"); axis square; title("specc K=2");
nexttile
scatter(embed2(:,1), embed2(:,2), 20, label3, 'filled');
xlabel("DIM1"); ylabel("DIM2"); axis square; title("specc K=3");
nexttile
scatter(embed2(:,1), embed2(:,2), 20, label4, 'filled');
xlabel("DIM1"); ylabel("DIM2"); axis square; title("specc K=4");