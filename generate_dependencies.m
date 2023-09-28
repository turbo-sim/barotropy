%% Save dependencies of the create_barotropic_polynomials.m script
% This script is used because I keep all many MATLAB functions in a local
% folder and I want to make a copy of them in a single place so that others
% can have access to them when downloading the script from a git repo

% Initialize script
clear all
close all
clc

% Create folder to save dependencies
dependencies_folder = '/dependencies';
if not(isfolder(dependencies_folder))
    mkdir(dependencies_folder)
end

% Get a list of all the dependencies
script_name = fullfile(fileparts(cd), 'barotropic_model.m');
function_list = matlab.codetools.requiredFilesAndProducts(script_name);
function_list(strcmp(function_list, script_name)) = [];

% Make a copy of each function
for i = 1:numel(function_list)
    copyfile(function_list{i}, dependencies_folder)
end