close all
clear;
clc
% load('FigureBlayeredModel3DTraScheme110.mat')
load('FigureClayeredModel3D110.mat')

plotimage(squeeze(Vx(:,80,:))')
xlabel('x/dx')
ylabel('z/dz')
title(' ')

plotimage(squeeze(Vy(:,80,:))')
xlabel('x/dx')
ylabel('z/dz')
title(' ')

plotimage(squeeze(Vz(:,80,:))')
xlabel('x/dx')
ylabel('z/dz')
title(' ')

plotimage(squeeze(DeltaH(:,80,:))')
xlabel('x/dx')
ylabel('z/dz')
title(' ')

plotimage(squeeze(DeltaV(:,80,:))')
xlabel('x/dx')
ylabel('z/dz')
title(' ')