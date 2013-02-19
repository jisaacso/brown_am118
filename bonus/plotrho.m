function plotrho
clear all; close all; hold on;
for(i=0:.001:2)
    plot(i,rho(1,i,[0,1,2],[.5,1.5]))
    plot(i,rho(2,i,[0,1,2],[.5,1.5]))
    plot(i,rho(3,i,[0,1,2],[.5,1.5]))
end
end