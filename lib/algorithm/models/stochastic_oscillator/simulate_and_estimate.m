function [MSE,C,x,f] = simulate_and_estimate(dts,Tend)
%% simulate and estimate - Estimate results for one draw
%
% Syntax:
% [MSE,C] = simulate and estimate(dts,Tend)
%
% In:
% dts - Time discrteizations (in seconds)
% Tend - End time (in seconds)
%
% Out:
% MSE - Mean squared error for each dt
% C - Standard deviation estimate from the filtering estimate
%
% Description:
% A realization of a frequency trajectory in the specific band  is
% simulated and then the function 'simulate periodic data' is called
% to simulate realizations of stochastic oscillatory signals with
% given parameters. Noise is added to simulate noisy observations. The MSE
% for each dt/TR is captured and returned.
%

% Parameters
N = Tend/dts(1);
T = dts(1)*(0:N-1);
Qc = 0.51;
x0 = [0;1];

% Set up frequency trajectory
f = 0.1*cumsum(randn(1,N)); % random walk
f = 1./(1+exp(f)); % transformed
f = 0.5+1*f; % in range (Hz)

MSE=0; C=0;

% Simulate full periodic data
x = simulate_periodic_data(N,dts(1),f,Qc,x0);

end