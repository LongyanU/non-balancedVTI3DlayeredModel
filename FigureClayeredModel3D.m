%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%                June 21       %%%%%%%%%%%%%%%%%
%%%%%%%%%       VTI 3D model           %%%%%%%%%%%%%%%%%
% 14:57 begin

% 时间已过 1569.009869 秒。
clear
clc %%%%%%%
close all
% Elapsed time is 6.198007 seconds.  May 2,2017
nt=420;    % number of time steps
eps=.6;     % stability
isnap=1;    % snapshot sampling

nx=250;
nz=250;
ny=250;

% vp=ones(nz,nx)*2000;
% yibuxiu=ones(nz,nx)*0.13;
% delta=ones(nz,nx)*0.13;


v=zeros(nz,nx,ny);
yibuxiu=zeros(nz,nx,ny);
delta=zeros(nz,nx,ny);
for i=1:nx
    for j=1:ny
        for k=1:120
            v(i,j,k)=2500;
            yibuxiu(i,j,k)=0.1;
            delta(i,j,k)=0.1;
        end
    end
end

for i=1:nx
    for j=1:ny
        for k=121:140
            v(i,j,k)=3000;
            yibuxiu(i,j,k)=0.2;
            delta(i,j,k)=0.05;
        end
    end
end

for i=1:nx
    for j=1:ny
        for k=141:160
            v(i,j,k)=3300;
            yibuxiu(i,j,k)=0.15;
            delta(i,j,k)=0.04;
        end
    end
end

for i=1:nx
    for j=1:ny
        for k=161:nz
            v(i,j,k)=2800;
            yibuxiu(i,j,k)=0.11;
            delta(i,j,k)=0.06;
        end
    end
end

vp=v;

dx=10;  %calculate space increment

x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=110;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian




seis_record=zeros(nt,nx);


% Source location
zs=100;
xs=nz/2;
ys=nz/2;

h=dx;

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics

tic

p=zeros([nx ny nz]); Vx=p; Vz=p;Vy=p;
DeltaH=p;
DeltaV=p;

% New coeff New FD scheme
coeff=[ 1.56201, -0.293775, 0.101451, -0.0423455, 0.0185804, -0.00793351, 0.00308954, -0.00101341, 0.000244924, -0.0000318866]


for it=1:320,
    it
    %DeltaH/x
    DeltaHx=(DeltaH-circshift(DeltaH,[1 0 0]));
    
    %DeltaH/x
    DeltaHy=(DeltaH-circshift(DeltaH,[0 1 0]));
    
    %DeltaV/z
    DeltaVz=(DeltaV-circshift(DeltaV,[0 0 1]));
    
    
    Vx=Vx+dt*(DeltaHx)/h;
    Vy=Vy+dt*(DeltaHy)/h;
    Vz=Vz+dt*(DeltaVz)/h;
%     [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    Vxx=coeff(1)*(circshift(Vx,[ -1 0 0])-circshift(Vx,[ 0 0 0]))+...
        coeff(2)*(circshift(Vx,[ -2 0 0])-circshift(Vx,[ 1 0 0]))+...
        coeff(3)*(circshift(Vx,[ -3 0 0])-circshift(Vx,[ 2 0 0]))+...
        coeff(4)*(circshift(Vx,[ -4 0 0])-circshift(Vx,[ 3 0 0]))+...
        coeff(5)*(circshift(Vx,[ -5 0 0])-circshift(Vx,[ 4 0 0]))+...
        coeff(6)*(circshift(Vx,[ -6 0 0])-circshift(Vx,[ 5 0 0]))+...
        coeff(7)*(circshift(Vx,[ -7 0 0])-circshift(Vx,[ 6 0 0]))+...
        coeff(8)*(circshift(Vx,[ -8 0 0])-circshift(Vx,[ 7 0 0]))+...
        coeff(9)*(circshift(Vx,[ -9 0 0])-circshift(Vx,[ 8 0 0]))+...
        coeff(10)*(circshift(Vx,[ -10 0 0])-circshift(Vx,[ 9 0 0]));
    
    Vyy=coeff(1)*(circshift(Vy,[0 -1 0])-circshift(Vy,[0 0 0]))+...
        coeff(2)*(circshift(Vy,[0 -2 0])-circshift(Vy,[0 1 0]))+...
        coeff(3)*(circshift(Vy,[0 -3 0])-circshift(Vy,[0 2 0]))+...
        coeff(4)*(circshift(Vy,[0 -4 0])-circshift(Vy,[0 3 0]))+...
        coeff(5)*(circshift(Vy,[0 -5 0])-circshift(Vy,[0 4 0]))+...
        coeff(6)*(circshift(Vy,[0 -6 0])-circshift(Vy,[0 5 0]))+...
        coeff(7)*(circshift(Vy,[0 -7 0])-circshift(Vy,[0 6 0]))+...
        coeff(8)*(circshift(Vy,[0 -8 0])-circshift(Vy,[0 7 0]))+...
        coeff(9)*(circshift(Vy,[0 -9 0])-circshift(Vy,[0 8 0]))+...
        coeff(10)*(circshift(Vy,[0 -10 0])-circshift(Vy,[0 9 0]));
    
    Vzz=coeff(1)*(circshift(Vz,[0 0 -1])-circshift(Vz,[0 0 0]))+...
        coeff(2)*(circshift(Vz,[0 0 -2])-circshift(Vz,[0 0 1]))+...
        coeff(3)*(circshift(Vz,[0 0 -3])-circshift(Vz,[0 0 2]))+...
        coeff(4)*(circshift(Vz,[0 0 -4])-circshift(Vz,[0 0 3]))+...
        coeff(5)*(circshift(Vz,[0 0 -5])-circshift(Vz,[0 0 4]))+...
        coeff(6)*(circshift(Vz,[0 0 -6])-circshift(Vz,[0 0 5]))+...
        coeff(7)*(circshift(Vz,[0 0 -7])-circshift(Vz,[0 0 6]))+...
        coeff(8)*(circshift(Vz,[0 0 -8])-circshift(Vz,[0 0 7]))+...
        coeff(9)*(circshift(Vz,[0 0 -9])-circshift(Vz,[0 0 8]))+...
        coeff(10)*(circshift(Vz,[0 0 -10])-circshift(Vz,[0 0 9]));
    
    DeltaH=DeltaH+dt*vp.^2.*( (1+2*yibuxiu).*(Vxx+Vyy)+sqrt(1+2*delta).*Vzz)/h;
    DeltaV=DeltaV+dt*vp.^2.*(sqrt(1+2*delta).*(Vxx+Vyy)+Vzz)/h;
    
    DeltaH(125,125,100)=DeltaH(125,125,100)+src(it);
    DeltaV(125,125,100)=DeltaV(125,125,100)+src(it);
    
    if rem(it,isnap)== 0,
        imagesc(x,z,squeeze(Vx(110,:,:)) ) , axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(max(Vx)))))
        drawnow
    end
    
end

toc
%
save('FigureClayeredModel3D110.mat')