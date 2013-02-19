function [z] = hw1_2ODE(u,t)
% We are trying to solve:
% integal(sin(t)*exp(-t),t,4,8)
% 
% Sample commandline call:
% hw1_1(4,8,0,1/128,@hw1_1ODE)

z = sin(t)*exp(-t);