clear all;
clc;
currentFolder = pwd;
addpath(genpath(currentFolder))

% parameter setting
fobj = @Rosenbrock; % Objective function
Np = 30 ; % population size  (Np is set to an even number greater than 2)
N = 10 ; % number of chaotic samples 
Dim = 10 ; % Dimensions of the optimization problem
Varmin = -30*ones(1,Dim); % Lower bound of optimization problem
Varmax = 30*ones(1,Dim); % Upper bound of optimization problem
MaxFES = Dim*1e4 ; % Maximum number of evaluation functions
tic
[Best,fBest,history] = CEO(fobj,Np,Dim,Varmin,Varmax,N,MaxFES);
toc
semilogy(history)


