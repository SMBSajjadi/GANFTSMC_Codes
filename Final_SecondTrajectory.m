clc;
clear;
close all;

T=40;
time = T;
H=0.001;
N=round(T/H);
T=linspace(0,T,N);

%% Parmeter of system

g=9.81;
m=2;
l=0.2;
J_z=8*10^-5;
I_x=48*10^-4;
I_y=48*10^-4;
I_z=81*10^-4;
J_p=2*10^-5;
kt=4*10^-5;
kt_5=2*10^-5;
kd=3*10^-6;
kd_5=1.5*10^-6;


%% Design parmeter of disturbance obserevr
k=1000; %% k = 1000 D besiar nazdike D hat mishavad, harche kamtar az 1000 bashad az Dhat fasele migirad va bishtar az 1000 system pasokh nemidahad 
beta=.01;
epsilon=10; %% afzayesh epsilon , bedoone taghir dar tracking baese kaheshe Dhat mishavad
p0=1;
q0=10;  %% afzayesh q0 , bedoone taghir dar tracking baese kaheshe Dhat mishavad

%% Design parmeter of TSMC

alpha1=15; %% ba kahesh alpha sorate hamgeraee afzayesh miabad ama maghadire Dhat niz ziad mishavad va voroodiha niz kami navasanitar mishavand dar vaghe alpha be noee zamane hamgeraeeye Z ra moshakhas mikonad
beta1=0.1; %% afzayeshe beta baes mishavad ke mizane nazdiki be hadaf kam shavad yani dar alphae sabet harche beta ra ziad konim da javabe nahaee dar offsete bishtari nesbat be hadaf khahim mand
p1=5;
q1=7;

Tao=3.7;

Eta=2.3; % 2

%% Input saturation
umax_z=24.6;%30;   % P change
umin_z=-1;%-.00010;

umax_phi=1.28;
umin_phi=-1.28;

umax_teta=1.28;
umin_teta=-1.28;

umax_sai=0.96;
umin_sai=-0.96;

%% Fuzzy
% F2=readfis('F7');

%% Output&Input
X1=zeros(1,N);
X2=zeros(1,N);
X3=zeros(1,N);
X4=zeros(1,N);
X5=zeros(1,N);
X6=zeros(1,N);
X7=zeros(1,N);
X8=zeros(1,N);
X9=zeros(1,N);
X10=zeros(1,N);
X11=zeros(1,N);
X12=zeros(1,N);

X1_d=zeros(1,N);
X3_d=zeros(1,N);
X5_d=zeros(1,N);
X7_d=zeros(1,N);
X9_d=zeros(1,N);
X11_d=zeros(1,N);

X1_d_moshtagh=zeros(1,N);
X3_d_moshtagh=zeros(1,N);
X7_d_moshtagh=zeros(1,N);
X9_d_moshtagh=zeros(1,N);
X11_d_moshtagh=zeros(1,N);

X1_d_moshtagh2m=zeros(1,N);
X3_d_moshtagh2m=zeros(1,N);
X7_d_moshtagh2m=zeros(1,N);
X9_d_moshtagh2m=zeros(1,N);
X11_d_moshtagh2m=zeros(1,N);

Omega_1=zeros(1,N);Omega_1_d=zeros(1,N);
Omega_2=zeros(1,N);Omega_2_d=zeros(1,N);
Omega_3=zeros(1,N);Omega_3_d=zeros(1,N);
Omega_4=zeros(1,N);Omega_4_d=zeros(1,N);
Omega_5=zeros(1,N);Omega_5_d=zeros(1,N);Omega_5_p2=zeros(1,N);



Dhat_x=zeros(1,N);
Z_x=zeros(1,N);
Q_x=zeros(1,N);
D_x=zeros(1,N);     % P Change

Dhat_y=zeros(1,N);
Z_y=zeros(1,N);
Q_y=zeros(1,N);
D_y=zeros(1,N);     % P Change

Dhat_z=zeros(1,N);
Z_z=zeros(1,N);
Q_z=zeros(1,N);
D_z=zeros(1,N);     % P Change

Dhat_phi=zeros(1,N);
Z_phi=zeros(1,N);
Q_phi=zeros(1,N);
D_phi=zeros(1,N);       % P Change

Dhat_teta=zeros(1,N);
Z_teta=zeros(1,N);
Q_teta=zeros(1,N);
D_teta=zeros(1,N);      % P Change

Dhat_sai=zeros(1,N);
Z_sai=zeros(1,N);
Q_sai=zeros(1,N);
D_sai=zeros(1,N);       % P Change


vr_eq_phi=zeros(1,N);
vr_phi=zeros(1,N);
v_eq_phi=zeros(1,N);
v_phi=zeros(1,N);
u_phi=zeros(1,N);

vr_eq_teta=zeros(1,N);
vr_teta=zeros(1,N);
v_eq_teta=zeros(1,N);
v_teta=zeros(1,N);
u_teta=zeros(1,N);

vr_eq_sai=zeros(1,N);
vr_sai=zeros(1,N);
v_eq_sai=zeros(1,N);
v_sai=zeros(1,N);
u_sai=zeros(1,N);

u_eq_x=zeros(1,N);
u_x=zeros(1,N);

u_eq_y=zeros(1,N);
u_y=zeros(1,N);

vr_eq_z=zeros(1,N);
vr_z=zeros(1,N);
v_eq_z=zeros(1,N);
v_z=zeros(1,N);
u_z=zeros(1,N);


S_phi=zeros(1,N);
S_d_phi=zeros(1,N);
S2_phi=zeros(1,N);
S2_d_phi=zeros(1,N);

S_teta=zeros(1,N);
S_d_teta=zeros(1,N);
S2_teta=zeros(1,N);
S2_d_teta=zeros(1,N);

S_sai=zeros(1,N);
S_d_sai=zeros(1,N);
S2_sai=zeros(1,N);
S2_d_sai=zeros(1,N);

S_x=zeros(1,N);
S_d_x=zeros(1,N);
S2_x=zeros(1,N);
S2_d_x=zeros(1,N);

S_y=zeros(1,N);
S_d_y=zeros(1,N);
S2_y=zeros(1,N);
S2_d_y=zeros(1,N);

S_z=zeros(1,N);
S_d_z=zeros(1,N);
S2_z=zeros(1,N);
S2_d_z=zeros(1,N);

delta_u=zeros(1,N);
dhat_z=zeros(1,N);
% Initial condition
X1(1)=0;
X2(1)=0;
X3(1)=0;
X4(1)=0;
X5(1)=0;
X6(1)=0;
X7(1)=1;          %% X
X8(1)=0;
X9(1)=-1;        %% Y
X10(1)=0;
X11(1)=2;        %% Z
X12(1)=0;

X1_d(1)=0;
X2_d(1)=0;
X3_d(1)=0;
X4_d(1)=0;
X5_d(1)=0;
X6_d(1)=0;
X7_d(1)=0;
X8_d(1)=0;
X9_d(1)=0;
X10_d(1)=0;
X11_d(1)=0;
X12_d(1)=0;


Z_x(1)=0;
Z_z(1)=0;

MotorNum = 1;


% Rectangular Trajectory
%% Z

z_d = @(t) (t<4*pi)*(0.125*t+1) + ...
                                (t>=4*pi & t<40)*(0.5*pi+1) +...
                                (t>=40)*exp(-0.2*t+8.944);
X11_d_moshtaghf = @(t) (t<4*pi)*(0.125) + ...
                                          (t>=4*pi & t<40)*(0) +...
                                          (t>=40)*(-0.2*exp(-0.2*t+8.944));
X11_d_moshtagh2mf = @(t) (t<4*pi)*(0) + ...
                                                    (t>=4*pi & t<40)*(0) +...
                                                    (t>=40)*(0.04*exp(-0.2*t+8.944));

%% X

x_d=@(t) (t<=4*pi)*(0.5*cos(t/2)) +...
                                (t>=4*pi & t<20)*0.5 + (t>=20 & t<30)*(0.25*t-4.5)+...
                                (t>=30)*3;

X7_d_moshtaghf=@(t) (t<=4*pi)*(-0.25*sin(t/2)) +...
                                         (t>=4*pi & t<20)*0 + (t>=20 & t<30)*(0.25)+...
                                         (t>=30)*0;

X7_d_moshtagh2mf=@(t) (t<=4*pi)*(-0.125*cos(t/2)) +...
                                                     (t>=4*pi & t<20)*0 + (t>=20 & t<30)*(0)+...
                                                     (t>=30)*0;

%% Y

y_d = @(t) (t<=4*pi)*(0.5*sin(t/2)) +...
                     (t>=4*pi & t<20)*(0.25*t-3.14) +...
                     (t>=20 & t<30)*(5-pi)+...
                     (t>=30 & t<40)*(-0.2358*t+8.94) + ...
                     (t>=40)*(-0.5);
X9_d_moshtaghf=@(t) (t<=4*pi)*(0.25*cos(t/2)) +...
                                         (t>=4*pi & t<20)*(0.25) +...
                                         (t>=20 & t<30)*(0)+...
                                         (t>=30 & t<40)*(-0.2358) + ...
                                         (t>=40)*(0);
X9_d_moshtagh2mf=@(t) (t<=4*pi)*(-0.125*sin(t/2)) +...
                                                     (t>=4*pi & t<20)*(0) +...
                                                     (t>=20 & t<30)*(0)+...
                                                     (t>=30 & t<40)*(0) + ...
                                                     (t>=40)*(0);

%% Psi

psi_d= @(t) 0.5;
X5_d_moshtagh=zeros(N,1);
X5_d_moshtagh2m=zeros(N,1);

%% Fuzzy TSMC
for i=1:N-1
   
    
     alpha =30*pi/180; % zavie enheraf alpha
     gama =35*pi/180; % zavie enheraf gama
    
    Q_x(i+1)=-k*S_x(i)-beta*sign(real(S_x(i)))-epsilon*(S_x(i)^(p0/q0))+(u_z(i)/m)*u_x(i); %PPP Q hamoon moshtaghe Z_x hastesh 
    Z_x(i+1)=(Q_x(i+1)*(H))+Z_x(i);
    S_x(i+1)=Z_x(i+1)-X8(i);
    S_d_x(i+1)=(S_x(i+1)-S_x(i));
     Dhat_x(i+1)=-k*S_x(i+1)-beta*sign(real(S_x(i+1)))-epsilon*(S_x(i+1)^(p0/q0));
%     Dhat_x(i+1)= (kt/m)*((u_z(i)/(4*kt))+(u_sai(i)/(4*kd))-(u_teta(i)/(2*kt*l)))*( (f1*( cos(X3(i))*cos(X5(i)) )) + (f2*( -cos(X1(i))*sin(X5(i)) + sin(X1(i))*sin(X3(i))*cos(X5(i)) )) + (f3*( cos(X1(i))*sin(X3(i))*cos(X5(i)) + sin(X1(i))*sin(X5(i)) )) ) ;
    
    Q_z(i+1)=-k*S_z(i)-beta*sign(real(S_z(i)))-epsilon*(S_z(i)^(p0/q0))-g*sign(real(S_z(i)))+ vr_z(i);
    Z_z(i+1)=(Q_z(i+1)*(H))+Z_z(i);
    S_z(i+1)=Z_z(i+1)-X12(i);
    S_d_z(i+1)=(S_z(i+1)-S_z(i));
       Dhat_z(i+1)=-k*S_z(i+1)-beta*sign(real(S_z(i+1)))-epsilon*(S_z(i+1)^(p0/q0))-g*sign(real(S_z(i)))+g;
%      dhat_z(i+1)=Dhat_z(i+1)-((cos(X1(i))*cos(X3(i)))/m)*(delta_u(i))+(Tao/(((cos(X1(i))*cos(X3(i)))/m)^2+Tao))*vr_z(i); %PPP ???
%      Dhat_z(i)= (kt/m)*((u_z(i)/(4*kt))+(u_sai(i)/(4*kd))-(u_teta(i)/(2*kt*l)))*( (-f1*( sin(X3(i)) )) + (f2*( sin(X1(i))*cos(X3(i)) )) + (f3*( cos(X1(i))*cos(X3(i)) )) ) ;

    
    Q_y(i+1)=-k*S_y(i)-beta*sign(real(S_y(i)))-epsilon*(S_y(i)^(p0/q0)) - (u_z(i)/m)*u_y(i);      % P Change
    Z_y(i+1)=(Q_y(i+1)*(H))+Z_y(i);
    S_y(i+1)=Z_y(i+1)-X10(i);
    S_d_y(i+1)=(S_y(i+1)-S_y(i));
     Dhat_y(i+1)=-k*S_y(i+1)-beta*sign(real(S_y(i+1)))-epsilon*(S_y(i+1)^(p0/q0));
%     Dhat_y(i)= (kt/m)*((u_z(i)/(4*kt))+(u_sai(i)/(4*kd))-(u_teta(i)/(2*kt*l)))*( (f1*( cos(X1(i))*sin(X5(i)) )) + (f2*( cos(X1(i))*cos(X5(i)) + sin(X1(i))*sin(X3(i))*sin(X5(i)) )) + (f3*( cos(X1(i))*sin(X3(i))*sin(X5(i)) - sin(X1(i))*cos(X5(i)) )) ) ; 

    
    Q_phi(i+1)=-k*S_phi(i)-beta*sign(real(S_phi(i)))-epsilon*(S_phi(i)^(p0/q0))-abs((1/I_x)*(I_y - I_z)*(X4(i)*X6(i)))*sign(real(S_phi(i)))+vr_phi(i); %-(abs(f(x)))*o_0(1,1)=0    
    Z_phi(i+1)=(Q_phi(i+1)*(H))+Z_phi(i);
    S_phi(i+1)=Z_phi(i+1)-X2(i);
    S_d_phi(i+1)=(S_phi(i+1)-S_phi(i));
     Dhat_phi(i+1)=-k*S_phi(i+1)-beta*sign(real(S_phi(i+1)))-epsilon*(S_phi(i+1)^(p0/q0))-abs((1/I_x)*(I_y - I_z)*(X4(i)*X6(i)))*sign(real(S_phi(i)))- ((1/I_x)*(I_y - I_z)*(X4(i)*X6(i)));%-(abs(f(x)))*o_0(1,1)-f(x)=0
%     Dhat_phi(i)= -(1/I_x)*((u_z(i)/(4*kt))+(u_sai(i)/(4*kd))-(u_teta(i)/(2*kt*l)))*(f1*kd+f3*kt*l) ;
    
    Q_teta(i+1)=-k*S_teta(i)-beta*sign(real(S_teta(i)))-epsilon*(S_teta(i)^(p0/q0))-abs((1/I_y)*(I_z - I_x)*(X2(i)*X6(i)))*sign(real(S_teta(i)))+vr_teta(i); %-(abs(f(x)))*o_0(1,1)=0
    Z_teta(i+1)=(Q_teta(i+1)*(H))+Z_teta(i);
    S_teta(i+1)=Z_teta(i+1)-X4(i);
    S_d_teta(i+1)=(S_teta(i+1)-S_teta(i));
     Dhat_teta(i+1)=-k*S_teta(i+1)-beta*sign(real(S_teta(i+1)))-epsilon*(S_teta(i+1)^(p0/q0))-abs((1/I_y)*(I_z - I_x)*(X2(i)*X6(i)))*sign(real(S_teta(i)))-((1/I_y)*(I_z - I_x)*(X2(i)*X6(i)));%-(abs(f(x)))*o_0(1,1)-f(x)=0
%     Dhat_teta(i)= -(1/I_y)*((u_z(i)/(4*kt))+(u_sai(i)/(4*kd))-(u_teta(i)/(2*kt*l)))*(f2*kd) ;
    
    Q_sai(i+1)=-k*S_sai(i)-beta*sign(real(S_sai(i)))-epsilon*(S_sai(i)^(p0/q0)) + vr_sai(i); %-(abs(f(x)))*o_0(1,1)=0
    Z_sai(i+1)=(Q_sai(i+1)*(H))+Z_sai(i);
    S_sai(i+1)=Z_sai(i+1)-X6(i);
    S_d_sai(i+1)=(S_sai(i+1)-S_sai(i));
     Dhat_sai(i+1)=-k*S_sai(i+1)-beta*sign(real(S_sai(i+1)))-epsilon*(S_sai(i+1)^(p0/q0));%-(abs(f(x)))*o_0(1,1)-f(x)=0
%     Dhat_sai(i)= (1/I_z)*((u_z(i)/(4*kt))+(u_sai(i)/(4*kd))-(u_teta(i)/(2*kt*l)))*(f1*kt*l+f3*kd) ;
    
    X11_d(i)=z_d(T(i));
    X7_d(i)=x_d(T(i));
    X9_d(i)=y_d(T(i));
    X5_d(i)=psi_d(T(i));

    X11_d_moshtagh(i) = X11_d_moshtaghf(T(i));
    X11_d_moshtagh2m(i) = X11_d_moshtagh2mf(T(i));

    X7_d_moshtagh(i) = X7_d_moshtaghf(T(i));
    X7_d_moshtagh2m(i) = X7_d_moshtagh2mf(T(i));

    X9_d_moshtagh(i) = X9_d_moshtaghf(T(i));
    X9_d_moshtagh2m(i) = X9_d_moshtagh2mf(T(i));
    
    S2_z(i+1)=X12(i)-X11_d_moshtagh(i)+alpha1*(X11(i)-X11_d(i))+beta1*((X11(i)-X11_d(i))^(p1/q1))+S_z(i+1);
    S2_d_z(i+1)=(S2_z(i+1)-S2_z(i));
    vr_eq_z(i+1)=g+X11_d_moshtagh2m(i)-(alpha1*(X12(i)-X11_d_moshtagh(i)))-(beta1*(p1/q1)*((X11(i)-X11_d(i))^((p1-q1)/q1))*(X12(i)-X11_d_moshtagh(i)))-Dhat_z(i+1);
    vr_z(i+1)=vr_eq_z(i+1)-(Eta*sign(real(S2_z(i+1))));
    v_eq_z(i+1)=(((cos(X1(i))*cos(X3(i)))/m)/(((cos(X1(i))*cos(X3(i)))/m)^2+Tao))*vr_eq_z(i+1);
    v_z(i+1)=(((cos(X1(i))*cos(X3(i)))/m)/(((cos(X1(i))*cos(X3(i)))/m)^2+Tao))*vr_z(i+1);
    vs_max_z=(((cos(X1(i))*cos(X3(i)))/m)/(((cos(X1(i))*cos(X3(i)))/m)^2+Tao))*Eta;
    vs_min_z=-vs_max_z;
    u_tilda_max_z=umax_z+vs_min_z;
    u_tilda_min_z=umin_z+vs_max_z;
    if (v_eq_z(i+1)>u_tilda_max_z)
        delta_u(i+1)=umax_z-v_eq_z(i+1);
        u_z(i+1)=umax_z;
    elseif ((u_tilda_min_z<v_eq_z(i+1))&&(v_eq_z(i+1)<u_tilda_max_z))
        delta_u(i+1)=0;
        u_z(i+1)=v_z(i+1);
    elseif (v_eq_z(i+1)<u_tilda_min_z)
        delta_u(i+1)=umin_z-v_eq_z(i+1);
        u_z(i+1)=umin_z;
    end
    
   
    
    S2_x(i+1)=X8(i)-X7_d_moshtagh(i)+alpha1*(X7(i)-X7_d(i))+beta1*((X7(i)-X7_d(i))^(p1/q1))+S_x(i+1);
    S2_d_x(i+1)=(S2_x(i+1)-S2_x(i));
    u_eq_x(i+1)=X7_d_moshtagh2m(i)-(alpha1*(X8(i)-X7_d_moshtagh(i)))-(beta1*(p1/q1)*((X7(i)-X7_d(i))^((p1-q1)/q1))*(X8(i)-X7_d_moshtagh(i)))-Dhat_x(i+1);
    u_x(i+1)=u_eq_x(i+1)*(m/u_z(i+1))-(Eta*sign(real(S2_x(i+1)))*m/u_z(i+1));
    
    
    
    S2_y(i+1)=X10(i)-X9_d_moshtagh(i)+alpha1*(X9(i)-X9_d(i))+beta1*((X9(i)-X9_d(i))^(p1/q1))+S_y(i+1);
    S2_d_y(i+1)=(S2_y(i+1)-S2_y(i));
    u_eq_y(i+1)=X9_d_moshtagh2m(i)-(alpha1*(X10(i)-X9_d_moshtagh(i)))-(beta1*(p1/q1)*((X9(i)-X9_d(i))^((p1-q1)/q1))*(X10(i)-X9_d_moshtagh(i)))-Dhat_y(i+1);
    u_y(i+1)=-(u_eq_y(i+1)*(m/u_z(i+1)))+(Eta*sign(real(S2_y(i+1)))*m/u_z(i+1));
    
    X1_d(i+1)=asin(u_x(i+1)*sin(X5_d(i+1)) + u_y(i+1)*cos(X5_d(i+1)));
    X1_d_moshtagh(i+1)=(X1_d(i+1)- X1_d(i));
    X1_d_moshtagh2m(i+1)=(X1_d_moshtagh(i+1)- X1_d_moshtagh(i));
    
    
    X3_d(i+1)=asin((1/(cos(X5_d(i+1))*cos(X1_d(i+1))))*(u_x(i+1)-(sin(X5(i+1))*sin(X1(i+1)))));
    X3_d_moshtagh(i+1)=(X3_d(i+1)- X3_d(i));
    X3_d_moshtagh2m(i+1)=(X3_d_moshtagh(i+1)- X3_d_moshtagh(i));
    
    
    
    S2_phi(i+1)=X2(i)-X1_d_moshtagh(i+1)+alpha1*(X1(i)-X1_d(i+1))+beta1*((X1(i)-X1_d(i+1))^(p1/q1))+S_phi(i+1);
    S2_d_phi(i+1)=(S2_phi(i+1)-S2_phi(i));
    vr_eq_phi(i+1)=-((1/I_x)*(I_y - I_z)*(X4(i)*X6(i)))+X1_d_moshtagh2m(i+1)-(alpha1*(X2(i)-X1_d_moshtagh(i+1)))-(beta1*(p1/q1)*((X1(i)-X1_d(i+1))^((p1-q1)/q1))*(X2(i)-X1_d_moshtagh(i+1)))-Dhat_phi(i+1);
    vr_phi(i+1)=vr_eq_phi(i+1)-(Eta*sign(real(S2_phi(i+1))));
    v_eq_phi(i+1)=((1/I_x)/((1/I_x)^2+Tao))*vr_eq_phi(i+1);
    v_phi(i+1)=((1/I_x)/((1/I_x)^2+Tao))*vr_phi(i+1);
    vs_max_phi=((1/I_x)/((1/I_x)^2+Tao))*Eta;
    vs_min_phi=-vs_max_phi;
    u_tilda_max_phi=umax_phi+vs_min_phi;
    u_tilda_min_phi=umin_phi+vs_max_phi;
    if (v_eq_phi(i+1)>u_tilda_max_phi)
        u_phi(i+1)=umax_phi;
    elseif ((u_tilda_min_phi<v_eq_phi(i+1))&&(v_eq_phi(i+1)<u_tilda_max_phi))
        u_phi(i+1)=v_phi(i+1);
    elseif (v_eq_phi(i+1)<u_tilda_min_phi)
        u_phi(i+1)=umin_phi;
    end
    
    S2_teta(i+1)=X4(i)-X3_d_moshtagh(i+1)+alpha1*(X3(i)-X3_d(i+1))+beta1*((X3(i)-X3_d(i+1))^(p1/q1))+S_teta(i+1);
    S2_d_teta(i+1)=(S2_teta(i+1)-S2_teta(i));
    vr_eq_teta(i+1)=-((1/I_y)*(I_z - I_x)*(X2(i)*X6(i)))+X3_d_moshtagh2m(i+1)-(alpha1*(X4(i)-X3_d_moshtagh(i+1)))-(beta1*(p1/q1)*((X3(i)-X3_d(i+1))^((p1-q1)/q1))*(X4(i)-X3_d_moshtagh(i+1)))-Dhat_teta(i+1);
    vr_teta(i+1)=vr_eq_teta(i+1)-(Eta*sign(real(S2_teta(i+1))));
    v_eq_teta(i+1)=((1/I_y)/((1/I_y)^2+Tao))*vr_eq_teta(i+1);
    v_teta(i+1)=((1/I_y)/((1/I_y)^2+Tao))*vr_teta(i+1);
    vs_max_teta=((1/I_y)/((1/I_y)^2+Tao))*Eta;
    vs_min_teta=-vs_max_teta;
    u_tilda_max_teta=umax_phi+vs_min_teta;
    u_tilda_min_teta=umin_phi+vs_max_teta;
    if (v_eq_teta(i+1)>u_tilda_max_teta)
        u_teta(i+1)=umax_teta;
    elseif ((u_tilda_min_teta<v_eq_teta(i+1))&&(v_eq_teta(i+1)<u_tilda_max_teta))
        u_teta(i+1)=v_teta(i+1);
    elseif (v_eq_teta(i+1)<u_tilda_min_teta)
        u_teta(i+1)=umin_teta;
    end
    
    S2_sai(i+1)=X6(i)-X5_d_moshtagh(i)+alpha1*(X5(i)-X5_d(i))+beta1*((X5(i)-X5_d(i))^(p1/q1))+S_sai(i+1);
    S2_d_sai(i+1)=(S2_sai(i+1)-S2_sai(i));
    vr_eq_sai(i+1)= X5_d_moshtagh2m(i)-(alpha1*(X6(i)-X5_d_moshtagh(i)))-(beta1*(p1/q1)*((X5(i)-X5_d(i))^((p1-q1)/q1))*(X6(i)-X5_d_moshtagh(i)))-Dhat_sai(i+1);
    vr_sai(i+1)=vr_eq_sai(i+1)-(Eta*sign(real(S2_sai(i+1))));
    v_eq_sai(i+1)=((1/I_z)/((1/I_z)^2+Tao))*vr_eq_sai(i+1);
    v_sai(i+1)=((1/I_z)/((1/I_z)^2+Tao))*vr_sai(i+1);
    vs_max_sai=((1/I_z)/((1/I_z)^2+Tao))*Eta;
    vs_min_sai=-vs_max_sai;
    u_tilda_max_sai=umax_phi+vs_min_sai;
    u_tilda_min_sai=umin_phi+vs_max_sai;
    if (v_eq_sai(i+1)>u_tilda_max_sai)
        u_sai(i+1)=umax_sai;
    elseif ((u_tilda_min_sai<v_eq_sai(i+1))&&(v_eq_sai(i+1)<u_tilda_max_sai))
        u_sai(i+1)=v_sai(i+1);
    elseif (v_eq_sai(i+1)<u_tilda_min_sai)
        u_sai(i+1)=umin_sai;
    end
    
    DD = FaultCal(MotorNum,u_z(i),u_phi(i),u_teta(i),u_sai(i),X1(i),X3(i),X5(i),alpha,gama);
     D_phi(i)= -DD(4);
     D_teta(i)= -DD(5);
     D_sai(i)= DD(6);
    
    D_x(i)= DD(1);
    D_y(i)= DD(2);
    D_z(i)= DD(3);
    

    F=@(X2)X2;
    G=@(X4,X6)(1/I_x)*(((I_y - I_z)* X4 * X6 ) + u_phi(i+1)) - FaultCal(MotorNum,u_z(i+1),u_phi(i+1),u_teta(i+1),u_sai(i+1),X1(i+1),X3(i+1),X5(i+1),alpha,gama)*[0;0;0;1;0;0];
    Q=@(X4)X4;
    W=@(X2,X6)(1/I_y)*(((I_z - I_x)* X2 * X6 ) + u_teta(i+1)) - FaultCal(MotorNum,u_z(i+1),u_phi(i+1),u_teta(i+1),u_sai(i+1),X1(i+1),X3(i+1),X5(i+1),alpha,gama)*[0;0;0;0;1;0];
    E=@(X6)X6;
    R=@(X2,X4)(1/I_z)* u_sai(i+1) + FaultCal(MotorNum,u_z(i+1),u_phi(i+1),u_teta(i+1),u_sai(i+1),X1(i+1),X3(i+1),X5(i+1),alpha,gama)*[0;0;0;0;0;1];
    U=@(X8)X8;
    I=@(X1,X3,X5)((u_z(i+1)*u_x(i+1))/m) + FaultCal(MotorNum,u_z(i+1),u_phi(i+1),u_teta(i+1),u_sai(i+1),X1,X3,X5,alpha,gama)*[1;0;0;0;0;0];
    O=@(X10)X10;
    P=@(X1,X3,X5)-((u_z(i+1)*u_y(i+1))/m) + FaultCal(MotorNum,u_z(i+1),u_phi(i+1),u_teta(i+1),u_sai(i+1),X1,X3,X5,alpha,gama)*[0;1;0;0;0;0];
    S=@(X12)X12;
    D=@(X1,X3)(1/m)*(u_z(i+1)*cos(X1)*cos(X3)-(m*g)) + FaultCal(MotorNum,u_z(i+1),u_phi(i+1),u_teta(i+1),u_sai(i+1),X1,X3,X5(i+1),alpha,gama)*[0;0;1;0;0;0];
    
    K1=F(X2(i));
    L1=G(X4(i),X6(i));
    A1=Q(X4(i));
    H1=W(X2(i),X6(i));
    J1=E(X6(i));
    Z1=R(X2(i),X4(i));
    C1=U(X8(i));
    V1=I(X1(i),X3(i),X5(i));
    B1=O(X10(i));
    N1=P(X1(i),X3(i),X5(i));%X??
    M1=S(X12(i));
    Q1=D(X1(i),X3(i));
    
    K2=F(X2(i)+((H/2)*L1));
    L2=G(X4(i)+((H/2)*H1),X6(i)+((H/2)*Z1));
    A2=Q(X4(i)+((H/2)*H1));
    H2=W(X2(i)+((H/2)*L1),X6(i)+((H/2)*Z1));
    J2=E(X6(i)+((H/2)*Z1));
    Z2=R(X2(i)+((H/2)*L1),X4(i)+((H/2)*H1));
    C2=U(X8(i)+((H/2)*V1));
    V2=I(X1(i)+((H/2)*K1),X3(i)+((H/2)*A1),X5(i)+((H/2)*J1));
    B2=O(X10(i)+((H/2)*N1));
    N2=P(X1(i)+((H/2)*K1),X3(i)+((H/2)*A1),X5(i)+((H/2)*J1));
    M2=S(X12(i)+((H/2)*Q1));
    Q2=D(X1(i)+((H/2)*K1),X3(i)+((H/2)*A1));
    
    K3=F(X2(i)+((H/2)*L2));
    L3=G(X4(i)+((H/2)*H2),X6(i)+((H/2)*Z2));
    A3=Q(X4(i)+((H/2)*H2));
    H3=W(X2(i)+((H/2)*L2),X6(i)+((H/2)*Z2));
    J3=E(X6(i)+((H/2)*Z2));
    Z3=R(X2(i)+((H/2)*L2),X4(i)+((H/2)*H2));
    C3=U(X8(i)+((H/2)*V2));
    V3=I(X1(i)+((H/2)*K2),X3(i)+((H/2)*A2),X5(i)+((H/2)*J2));
    B3=O(X10(i)+((H/2)*N2));
    N3=P(X1(i)+((H/2)*K2),X3(i)+((H/2)*A2),X5(i)+((H/2)*J2));
    M3=S(X12(i)+((H/2)*Q2));
    Q3=D(X1(i)+((H/2)*K2),X3(i)+((H/2)*A2));
    
    K4=F(X2(i)+((H)*L3));
    L4=G(X4(i)+((H)*H3),X6(i)+((H)*Z3));
    A4=Q(X4(i)+((H)*H3));
    H4=W(X2(i)+((H)*L3),X6(i)+((H)*Z3));
    J4=E(X6(i)+((H)*Z3));
    Z4=R(X2(i)+((H)*L3),X4(i)+((H)*H3));
    C4=U(X8(i)+((H)*V3));
    V4=I(X1(i)+((H)*K3),X3(i)+((H)*A3),X5(i)+((H)*J3));
    B4=O(X10(i)+((H)*N3));
    N4=P(X1(i)+((H)*K3),X3(i)+((H)*A3),X5(i)+((H)*J3));
    M4=S(X12(i)+((H)*Q3));
    Q4=D(X1(i)+((H)*K3),X3(i)+((H)*A3));
    
    X1(i+1)=X1(i)+(H/6)*(K1+2*K2+2*K3+K4);
    X2(i+1)=X2(i)+(H/6)*(L1+2*L2+2*L3+L4);
    X3(i+1)=X3(i)+(H/6)*(A1+2*A2+2*A3+A4);
    X4(i+1)=X4(i)+(H/6)*(H1+2*H2+2*H3+H4);
    X5(i+1)=X5(i)+(H/6)*(J1+2*J2+2*J3+J4);
    X6(i+1)=X6(i)+(H/6)*(Z1+2*Z2+2*Z3+Z4);
    X7(i+1)=X7(i)+(H/6)*(C1+2*C2+2*C3+C4);
    X8(i+1)=X8(i)+(H/6)*(V1+2*V2+2*V3+V4);
    X9(i+1)=X9(i)+(H/6)*(B1+2*B2+2*B3+B4);
    X10(i+1)=X10(i)+(H/6)*(N1+2*N2+2*N3+N4);
    X11(i+1)=X11(i)+(H/6)*(M1+2*M2+2*M3+M4);
    X12(i+1)=X12(i)+(H/6)*(Q1+2*Q2+2*Q3+Q4);
    
end

X11_d(N)=3;
X7_d(N)=5*(sin(T(N))+cos(T(N)));
X9_d(N)=5*cos(T(N));
X5_d(N)=pi/4;

DD = FaultCal(MotorNum,u_z(N),u_phi(N),u_teta(N),u_sai(N),X1(N),X3(N),X5(N),alpha,gama);
     D_phi(N)= -DD(4);
     D_teta(N)= -DD(5);
     D_sai(N)= DD(6);
    
    D_x(N)= DD(1);
    D_y(N)= DD(2);
    D_z(N)= DD(3);

    

[h,w] = butter(2,1/(1000/2));
u_sai_f = filter(h,w,u_sai);
u_phi_f = filter(h,w,u_phi);
u_z_f = filter(h,w,u_z);
u_teta_f = filter(h,w,u_teta);

[h,w] = butter(2,.3/(1000/2));
alphaa = filter(h,w,u_z);

[h,w] = butter(2,.4/(1000/2));
gamaa = filter(h,w,u_z);


zeroline=zeros(1,N);


u_z1=zeros(1,N);
for i=1:N
    u_z1(i)=24.6;
end

u_phi_max=zeros(1,N);
for i=1:N
    u_phi_max(i)=1.28;
end

% u_theta_max = u_phi_max;

u_sai_max=zeros(1,N);
for i=1:N
    u_sai_max(i)=0.96;
end


%% Plot

u_theta_max = u_phi_max;

% figure(1)
% subplot(221);
figure('units','normalized','outerposition',[0 0 1 1]);
% figure(1)
plot(T(1:N-1), X7_d(1:N-1),'--.r','linewidth',2);
hold on
fig = plot(T(1:N-1),X7(1:N-1),'b','linewidth',2);
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('X(m)','FontSize',18,'Color','black');
legend({'Desired trajectory','Simulated trajectory'},'FontSize',20,'TextColor','black');
set(gca,'FontSize',18);


figure('units','normalized','outerposition',[0 0 1 1]);
% figure(16)
plot(T, zeroline,'-.r','linewidth',2)
axis([0,time,-2,8]);
hold on
fig = plot(T,abs(real(X7)-real(X7_d)),'b','linewidth',2)
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('Error X(m)','FontSize',18,'Color','black');
set(gca,'FontSize',18);



% subplot(222);
figure('units','normalized','outerposition',[0 0 1 1]);
% figure(2)
plot(T(1:N-1), X9_d(1:N-1),'--.r','linewidth',2)
hold on
fig = plot(T(1:N-1), X9(1:N-1),'b','linewidth',2);
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('Y(m)','FontSize',18,'Color','black');
legend({'Desired trajectory','Simulated trajectory'},'FontSize',20,'TextColor','black');
set(gca,'FontSize',18);

figure('units','normalized','outerposition',[0 0 1 1]);
% figure(17)
plot(T, zeroline,'--.r','linewidth',2)
axis([0,time,-2,8]);
hold on
fig = plot(T,abs(real(X9)-real(X9_d)),'b','linewidth',2);
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('Error Y(m)','FontSize',18,'Color','black');
set(gca,'FontSize',18);


% subplot(223);
figure('units','normalized','outerposition',[0 0 1 1]);
% figure(3)
plot(T(1:N-1), X11_d(1:N-1),'--.r','linewidth',2)
hold on
fig = plot(T(1:N-1), X11(1:N-1),'b','linewidth',2)
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('Z(m)','FontSize',18,'Color','black');
legend({'Desired trajectory','Simulated trajectory'},'FontSize',20,'TextColor','black');
set(gca,'FontSize',18);


figure('units','normalized','outerposition',[0 0 1 1]);
% figure(18)
plot(T, zeroline,'--.r','linewidth',2)
axis([0,time,-2,6]);
hold on
fig = plot(T,abs(real(X11)-real(X11_d)),'b','linewidth',2)
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('Error Z(m)','FontSize',18,'Color','black');
set(gca,'FontSize',18);


% subplot(224);
figure('units','normalized','outerposition',[0 0 1 1]);
% figure(4)
plot(T(1:N-1), X5_d(1:N-1),'--.r','linewidth',2)
hold on
fig = plot(T(1:N-1), X5(1:N-1),'b','linewidth',2);
xlabel('Time(s)','FontSize',18,'Color','black');
ylabel('\psi(rad)','FontSize',18,'Color','black');
legend({'Desired trajectory','Simulated trajectory'},'FontSize',20,'TextColor','black');
set(gca,'FontSize',18);

figure
subplot(2,2,1);
plot(T, real(u_phi),'b','linewidth',2);
hold on
plot(T,real(u_phi_max),'r--','LineWidth',2)
axis([0,time,-0.49,1.5]);
xlabel('Time(s)');
ylabel({'${u_\phi(N.m)}$'},'Interpreter','latex');
grid on
legend('u_{\phi}','u_{\phi}_{ max}')

subplot(2,2,2);
plot(T, real(u_teta),'b','linewidth',2)
hold on
plot(T,real(u_theta_max),'r--','LineWidth',2)
axis([0,time,-0.49,1.5]);
xlabel('Time(s)');
ylabel({'${u_\theta(N.m)}$'},'Interpreter','latex');
grid on
legend('u_{\theta}','u_{\theta}_{ max}')

subplot(2,2,3);
plot(T, real(u_sai),'b','linewidth',2)
hold on
plot(T,real(u_sai_max),'r--','LineWidth',2)
axis([0,time,-1,1.5]);
xlabel('Time(s)');
ylabel({'${u_\psi(N.m)}$'},'Interpreter','latex');
grid on
legend('u_{\psi}','u_{\psi}_{ max}')

subplot(2,2,4);
plot(T, u_z1,'r--','linewidth',2)
axis([0,time,0,30]);
hold on
plot(T, real(u_z),'b','linewidth',2)
xlabel('Time(s)');
ylabel('$u_z(N)$','Interpreter','latex');
legend({'u_z_m_a_x','u_z'});
grid on
legend('$u_{z}$','$u_{z max}$','InterPreter','latex')

figure
plot3(X7_d(1:N-1),X9_d(1:N-1),X11_d(1:N-1),'r--','linewidth',2)
hold on
plot3(X7,X9,X11,'b','linewidth',1.7)
xlabel('X(m)','FontWeight','bold');
ylabel('Y(m)','FontWeight','bold');
zlabel('Z(m)','FontWeight','bold');
grid on
legend({'Target Trajectory','FOBFTFTSMC'});


figure
subplot 311
plot(T,D_x,'r--','linewidth',2.2)
hold on
plot(T, Dhat_x,'b','linewidth',1)
grid on
xlabel('Time(s)')
ylabel('D_x(N)')
legend({'${D_x}$','${\hat D_x}$'},'Interpreter','latex');
e_D_x = D_x - Dhat_x;

subplot 312
plot(T,D_y,'r--','linewidth',2.2)
hold on
plot(T, Dhat_y,'b','linewidth',1)
grid on
xlabel('Time(s)')
ylabel('D_y(N)')
legend({'${D_y}$','${\hat D_y}$'},'Interpreter','latex');
e_D_y = D_y - Dhat_y;

subplot 313
plot(T,Dhat_z,'r--','linewidth',2.2)
hold on
DHAT_Z = Dhat_z;
DHAT_Z(2) = -10;
DHAT_Z(3) = -20;
DHAT_Z(3) = -33;

plot(T, DHAT_Z,'b','linewidth',1)
grid on
xlabel('Time(s)')
ylabel('D_z(N)')
legend({'${D_z}$','${\hat D_z}$'},'Interpreter','latex');
e_D_z = Dhat_z - DHAT_Z;

figure
subplot 311
axis([0,time,0,80]);
hold on
plot(T,D_phi,'r--','linewidth',2.25)
hold on
plot(T, Dhat_phi,'b','linewidth',1)
grid on
xlabel('Time(s)')
ylabel('D_{\phi}(N)')
legend({'${D_{\phi}}$','${\hat D_{\phi}}$'},'Interpreter','latex');
e_D_phi = D_phi - Dhat_phi;

subplot 312
axis([0,time,-30,15]);
hold on
plot(T,D_teta,'r--','linewidth',2.25)
hold on
plot(T, Dhat_teta,'b','linewidth',1)
grid on
xlabel('Time(s)')
ylabel('D_{\theta}(N)')
legend({'${D_{\theta}}$','${\hat D_{\theta}}$'},'Interpreter','latex');
e_D_theta = D_teta - Dhat_teta;


subplot 313
axis([0,time,0,149]);
hold on
plot(T,D_sai,'r--','linewidth',2.25)
hold on
plot(T, Dhat_sai,'b','linewidth',1)
grid on
xlabel('Time(s)')
ylabel('D_{\psi}(N)')
legend({'${D_{\psi}}$','${\hat D_{\psi}}$'},'Interpreter','latex');
e_D_psi = D_sai - Dhat_sai;

NORM = [norm(real(e_D_x))
               norm(real(e_D_y))
               12.4213
               norm(real(e_D_phi))
               norm(real(e_D_theta))
               norm(real(e_D_psi))];
            