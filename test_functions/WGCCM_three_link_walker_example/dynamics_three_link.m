function [D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,dLfH]=dynamics_three_link(x,a)
% DYNAMICS_THREE_LINK    Model of three-link biped walker model.
%    [D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,DLFH] = DYNAMICS_THREE_LINK(X,
%      A) is the three-link
%    biped walking model. (x is of dimension 6)

% Eric Westervelt
% 20-Feb-2007 10:18:18

[r,m,Mh,Mt,L,g]=model_params_three_link;

[th3d,th1d,alpha,epsilon]=control_params_three_link;

th1=x(1); th2=x(2); th3=x(3);
dth1=x(4); dth2=x(5); dth3=x(6);

% D matrix
D=zeros(3);
D(1,1)=5/4*m*r^2+Mh*r^2+r^2*Mt;
D(1,2)=-1/2*r^2*cos(th1-th2)*m;
D(1,3)=Mt*r*L*cos(-th1+th3);
D(2,1)=-1/2*r^2*cos(th1-th2)*m;
D(2,2)=1/4*m*r^2;
D(3,1)=Mt*r*L*cos(-th1+th3);
D(3,3)=L^2*Mt;

% C matrix
C=zeros(3);
C(1,2)=-1/2*r^2*sin(th1-th2)*m*dth2;
C(1,3)=-Mt*r*L*sin(-th1+th3)*dth3;
C(2,1)=1/2*r^2*sin(th1-th2)*m*dth1;
C(3,1)=Mt*r*L*sin(-th1+th3)*dth1;

% G matrix
G=zeros(3,1);
G(1)=-3/2*r*sin(th1)*m*g-r*sin(th1)*Mh*g-r*sin(th1)*Mt*g;
G(2)=1/2*r*sin(th2)*m*g;
G(3)=-L*sin(th3)*Mt*g;

% B matrix
B=zeros(3,2);
B(1,1)=-1;
B(2,2)=-1;
B(3,1)=1;
B(3,2)=1;

% K matrix
K=zeros(2,4);
K(1,1)=156;
K(1,3)=25;
K(2,2)=110;
K(2,4)=21;

% dV matrix
dV=zeros(1,6);
dV(1,1)=184/77*th1+184/77*th2-10*dth1-10*dth2;
dV(1,2)=184/77*th1+184/77*th2-10*dth1-10*dth2;
dV(1,3)=4515147318722739/2251799813685248*th3-4515147318722739/4503599627370496*th3d-4515147318722739/4503599627370496*conj(th3d)-10*dth3;
dV(1,4)=-10*th1-10*th2+370/7*dth1+370/7*dth2;
dV(1,5)=-10*th1-10*th2+370/7*dth1+370/7*dth2;
dV(1,6)=-10*th3+5*th3d+4419157134357289/70368744177664*dth3+5*conj(th3d);

% dVl matrix
dVl=zeros(1,4);
dVl(1,1)=4515147318722739/2251799813685248*th3-4515147318722739/2251799813685248*conj(th3d)-10*dth3;
dVl(1,2)=184/77*th1+184/77*th2-10*dth1-10*dth2;
dVl(1,3)=-10*th3+10*conj(th3d)+4419157134357289/70368744177664*dth3;
dVl(1,4)=-10*th1-10*th2+370/7*dth1+370/7*dth2;

% Al matrix
Al=zeros(4,4);
Al(1,3)=1;
Al(2,4)=1;

% Bl matrix
Bl=zeros(4,2);
Bl(3,1)=1;
Bl(4,2)=1;
a01=a(1); a11=a(2); a21=a(3); a31=a(4);
a02=a(5); a12=a(6); a22=a(7); a32=a(8);

% Ha matrix
H=zeros(2,1);
H(1,1)=th3-a01-a11*th1-a21*th1^2-a31*th1^3;
H(2,1)=th2+th1-(a02+a12*th1+a22*th1^2+a32*th1^3)*(th1-th1d)*(th1+th1d);

% LfHa matrix
LfH=zeros(2,1);
LfH(1,1)=(-a11-2*a21*th1-3*a31*th1^2)*dth1+dth3;
LfH(2,1)=(1-(a12+2*a22*th1+3*a32*th1^2)*(th1-th1d)*(th1+th1d)-(a02+a12*th1+a22*th1^2+a32*th1^3)*(th1+th1d)-(a02+a12*th1+a22*th1^2+a32*th1^3)*(th1-th1d))*dth1+dth2;

% dLfHa matrix
dLfH=zeros(2,6);
dLfH(1,1)=(-2*a21-6*a31*th1)*dth1;
dLfH(1,4)=-a11-2*a21*th1-3*a31*th1^2;
dLfH(1,6)=1;
dLfH(2,1)=(-(2*a22+6*a32*th1)*(th1-th1d)*(th1+th1d)-2*(a12+2*a22*th1+3*a32*th1^2)*(th1+th1d)-2*(a12+2*a22*th1+3*a32*th1^2)*(th1-th1d)-2*a02-2*a12*th1-2*a22*th1^2-2*a32*th1^3)*dth1;
dLfH(2,4)=1-(a12+2*a22*th1+3*a32*th1^2)*(th1-th1d)*(th1+th1d)-(a02+a12*th1+a22*th1^2+a32*th1^3)*(th1+th1d)-(a02+a12*th1+a22*th1^2+a32*th1^3)*(th1-th1d);
dLfH(2,5)=1;
