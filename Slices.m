close all
clear;
clc
load('FigureBlayeredModel3DTraScheme110.mat')
aa=squeeze(Vx(100,80,:))';
bb=squeeze(Vy(100,80,:))';
cc=squeeze(Vz(100,80,:))';
% dd=squeeze(DeltaH(100,80,:))';
% ee=squeeze(DeltaV(100,80,:))';

load('FigureClayeredModel3D110.mat')
aaa=squeeze(Vx(100,80,:))';
bbb=squeeze(Vy(100,80,:))';
ccc=squeeze(Vz(100,80,:))';
dd=squeeze(DeltaH(100,80,:))';
% ee=squeeze(DeltaV(100,80,:))';

figure;plot(aa,'b')
hold on;plot(bb-1*10^-10,'b')
hold on;plot(cc-2*10^-10,'b')


hold on;plot(aaa,'k')
hold on;plot(bbb-1*10^-10,'k')
hold on;plot(ccc-2*10^-10,'k')

axis([1 250 -0.85*10^-9 0.85*10^-9])
legend('Traditional balanced SGFD scheme', 'Non-balanced SGFD scheme')

xlabel('z/dz')
ylabel('m/s')
