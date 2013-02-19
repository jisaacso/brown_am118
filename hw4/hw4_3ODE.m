function [z] = hw4_3ODE(u,t)

A = [-1 -1; 0 -1000];
z = A*u;