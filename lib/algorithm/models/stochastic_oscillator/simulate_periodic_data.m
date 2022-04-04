function [x] = simulate_periodic_data(N,dt,f,Qc,x0)
%% simulate and estimate ? Estimate results for one draw
%
% Syntax:
% simulate periodic data(N,dt,f,Qc,x0)

%
% In:
% N ? Number of steps to simulate
% dt ? Time discrteization step length (in seconds)
% f ? A vector of requencies in Hz
% Qc ? Spectral denisty of the dynamic noise term
% x0 ? Initial state
%
% Out:
% x -2xN-vector of simulated states
%
% Description:
% Simulate a stochastic oscillator with a given frequency trajectory.
% The model is set up as a continuous-time state space model, or a
% stochastic differential equation. 
%

%% Simulate a stochastic oscillator
% Allocate space for results
x = zeros(2,N);
% Initial state
if nargin < 5 || isempty(x0)
x(:,1) = randn(2,1);
else
x(:,1) = x0;
end
% The rest of the states
for k=2:N
% Dynamic model
F = [ 0 2*pi*f(k);
     -2*pi*f(k) 0];
L = [0;1];
% Discretize
[A,Q] = lti_disc(F,L,Qc,dt);
% Determine if stochastic oscillator or deterministic
    if Qc<eps
        x(:,k) = A*x(:,k-1);
    else
        x(:,k) = A*x(:,k-1) + chol(Q)'*randn(2,1);
    end
end

end