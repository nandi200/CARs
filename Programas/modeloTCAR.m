%Runge-Kutta para el modelo T-CAR
%dx/dt=A-a_0x-dxc_s C?LULAS NAIVE
%dy/dt=dxc_s-d_0y  C?LULAS INFECTADAS
%dc_s/dt=k-k_0c_s-k_1c_sz C?LULAS CANCEROSAS
%dz/dt=B+bc_s-b_0z-b_1zc_s C?LULAS INMUNES

%Constantes
A=1.5;
a1=0.01;
d=0.00001;
d2=0.003;
k=10;
k2=5;
k1=0.005;
B=2;
b=0.01;
b2=0.05;
b1=0.001;

%Define function handles
fX=@(t,X,CS)A-a1*X-d*X*CS;
fY=@(t,Y,X,CS)d*X*CS-d2*Y;
fCS=@(t,CS,Z)k-k2*CS-k1*CS*Z;
fZ=@(t,CS,Z)B+b*CS-b2*Z-b1*Z*CS;

%Initial Conditions
t(1)=0;
X(1)=149.7118;
CS(1)=0.9607;
Y(1)=1.9251;
Z(1)=38.8877;

%Step size
h=0.1;
tfinal=500;
N=ceil(tfinal/h);

%Update Loop
for i=1:N
    %Update time
    t(i+1)=t(i)+h;
    k1X=fX(t(i)      , X(i),       CS(i));
    k1Y=fY(t(i)      , X(i),        CS(i), Y(i)        );
    k1CS=fCS(t(i)    , Z(i),        CS(i));
    k1Z=fCS(t(i)     , Z(i),        CS(i));
    %k2
    k2X=fX(t(i)+h/2  , X(i)+h/2*k1X, CS(i)+h/2*k1CS);
    k2Y=fY(t(i)+h/2  , X(i)+h/2*k1X, CS(i)+h/2*k1CS, Y(i)+h/2*k1Y);
    k2CS=fCS(t(i)+h/2, Z(i)+h/2*k1Z, CS(i)+h/2*k1CS);
    k2Z=fZ(t(i)+h/2  , Z(i)+h/2*k1Z, CS(i)+h/2*k1CS);
    %k3
    k3X=fX(t(i)+h/2  , X(i)+h/2*k2X, CS(i)+h/2*k2CS);
    k3Y=fY(t(i)+h/2  , X(i)+h/2*k2X, CS(i)+h/2*k2CS, Y(i)+h/2*k2Y);
    k3CS=fCS(t(i)+h/2, Z(i)+h/2*k2Z, CS(i)+h/2*k2CS);
    k3Z=fZ(t(i)+h/2  , Z(i)+h/2*k2Z, CS(i)+h/2*k2CS);
    %k4
    k4X=fX(t(i)+h    , X(i)+h*k3X,  CS(i)+h*k3CS);
    k4Y=fY(t(i)+h    , X(i)+h*k3X,  CS(i)+h*k3CS    , Y(i)+h*k3Y);
    k4CS=fCS(t(i)+h  , Z(i)+h*k3Z,  CS(i)+h*k3CS);
    k4Z=fZ(t(i)+h    , Z(i)+h*k3Z,  CS(i)+h*k3CS);
    %x,y,cs,z
    X(i+1)=X(i)+h/6*(k1X+2*k2X+2*k3X+k4X);
    Y(i+1)=Y(i)+h/6*(k1Y+2*k2Y+2*k3Y+k4Y);
    CS(i+1)=CS(i)+h/6*(k1CS+2*k2CS+2*k3CS+k4CS);
    Z(i+1)=Z(i)+h/6*(k1Z+2*k2Z+2*k3Z+k4Z);
end

%Plot the solution
plot(t,CS)
hold on



    
