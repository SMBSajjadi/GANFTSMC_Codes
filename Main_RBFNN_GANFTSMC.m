clc;clear;close all;

%{
    By M. Sajjadi
%}

%% Main GA

% opt = optimoptions("ga","MaxGenerations",200,...
%                               "CrossoverFraction",0.2,...
%                               "MaxStallGenerations",50);
% 
% LB = 1e-3*ones(2,1);
% [Eta_K_Optimal,Feval,exitFlag] = ga(@(K)GeneticFunction(K),2,[],[],[],[],LB,[],[],opt);

%% FAULT ANGLES

T = pi/180;     % Conversion to Radians
ALPHA = 20*T;
GAMMA = 30*T;
BETA = 10*T;
FAULT_ANGLES = [ALPHA        BETA        GAMMA]';
m = 1;
g = 9.81;

%% Initialization

Ts = 0.001;   %% Must be 0.0001, at least, to be Fully Robust.
tMax = 50;
t = 0:Ts:tMax;
N = numel(t);
n = 6;  % Number of DOF (X)

x0 = [3 1 1 0 0 0 0 0 0 0 0 0];       %% Ics
x = zeros(2*n,N);                          %% Compelete State Vector
x(:,1) = x0;

s = @(x) sin(x);
c = @(x) cos(x);
S = 0.1*ones(n,N);

%% NFTSMC Parameters

%% Optimal Gains

% Eta =Eta_K_Optimal(1);
% K =Eta_K_Optimal(2);
Eta_K_Optimal = [0.512 8.242];

% Controller Parameters
b = [0.01 0.01 0.01 0.1 0.1 0.1];       
bPrime = [8 8 8 12 12 12];          
Landa = 2;
LandaPrime = 1.8;           



%% Non-Optimal Gains
% Eta = 2*[1 1 1 1 1 1];
% K = [15 15 15 15 15 15];
% b = [0.01 0.01 0.01 0.1 0.1 0.1];                    
% bPrime = [8 8 2 12 12 12];           
% Landa = 2;
% LandaPrime = 1.8;                

%% RBF Parameters

gammaX = 0.01;
gammaY = 0.01;
gammaZ = 0.01;

gammaPhi = 0.01;
gammaTheta = 0.01;
gammaPsi = 0.01;
GammaRBF = [gammaX, gammaY,gammaZ,gammaPhi,gammaTheta,gammaPsi]';

nKernel = 5;           % Number of Kernel Functions
BIJ = diag([0.1 0.1 0.1 0.4 0.4 0.4]);
CIJ = diag([0.01 0.01 0.01 0.001 0.001 0.001]);

bij = zeros(n,nKernel);
cij = bij;

for i=1:n
    
    bij(i,:) = BIJ(i,i);
    cij(i,:) = CIJ(i,i);
    
end

%% Reference Signals

CASE = 1;           %% Linear Trajectory
% CASE = 3;       %% BiLinear Trajectory
[XD, XDotD, XDoubleDotD] = setDesiredTrajectory(t,CASE,n);

phid = XD(4,:);
thetad = XD(5,:);
sayd = XD(6,:);

%% Control Signal Initialization

nU = 6;             %% iN cOnjunction with the Virtual Signals
u = zeros(nU,N);
u(3,1) = 0.1*m*g;
W = 100*ones(4,N);
W(:,1) = [120 100 120 100];
Error = zeros(n,N);
ErrorDot = 0.01*ones(size(Error));

%% Fault Estimation Initialization

FaultHat = 0.01*ones(n,N);
wHatX = 0.01*ones(nKernel,N);
wHatY = 0.01*ones(nKernel,N);
wHatZ= 0.01*ones(nKernel,N);
wHatPhi= 0.01*ones(nKernel,N);
wHatTheta= 0.01*ones(nKernel,N);
wHatPsi = 0.01*ones(nKernel,N);

%% Estimation of Fault Angles

nAngles = 3;            % Number of Faulty Angles
EstimatedFaultAngles = zeros(nAngles,N);

%% Fault Vectors

SystemStructuralFault = zeros(n,N);
Fi_Main_Fault = SystemStructuralFault;

%% Main Loop

tic;

for i=2:N  
    %% State Calculation
    
    [K1RK,W(:,i), SystemStructuralFault(:,i),wStar,Fi_Main_Fault(:,i)] =...
                                                                    Rotor2_Dynamic(t(i-1),...
                                                                                            x(:,i-1),...
                                                                                            u(:,i-1),...
                                                                                            FAULT_ANGLES);

    x(:,i) = stateCalculation(K1RK,x(:,i-1),u(:,i-1),Ts,t(i-1),FAULT_ANGLES);
    
    %% Virtual Control Design
    
    ux = u(1,i-1);
    uy = u(2,i-1);
    uz = u(3,i-1);
    
    thetad(i) = atan((ux*c(sayd(i))+uy*s(sayd(i)))/(g+uz));
    phid(i) = atan(c(thetad(i))*((ux*s(sayd(i))-uy*c(sayd(i)))/(g+uz)));
    
    %% Desired Roll/Pittch Angles
    
    XD(4:5,i) = [phid(i), thetad(i)]'; 

    %% RBFNNNTSMC Design
  
    [u(:,i),S(:,i),Error(:,i),ErrorDot(:,i)] = ...
        RBFNNTSMC_GENETIC(wStar,x(:,i),n,XD(:,i),XDotD(:,i),XDoubleDotD(:,i),...
        b, bPrime, Landa,LandaPrime, Eta_K_Optimal, FaultHat(:,i-1));

    %% RBF Estimator
    
    [wHatX(:,i),wHatY(:,i),wHatZ(:,i),wHatPhi(:,i),wHatTheta(:,i),wHatPsi(:,i), FaultHat(:,i)] = ...
    RBFNN_StructuralFaultEstimator(t,Error(:,i),ErrorDot(:,i),S(:,i),bij,cij,GammaRBF,...
    nKernel,wHatX(:,i-1),wHatY(:,i-1),wHatZ(:,i-1),...
    wHatPhi(:,i-1),wHatTheta(:,i-1),wHatPsi(:,i-1),bPrime,LandaPrime);
    
end

toc;

%% PlotResults

PlotComparedResults_GANFTSMC_NFTSMC()    
% plotResults(t,x,XD,XDotD,u,n,W,S)
% PlotEstimationResultsOfFaults(t,SystemStructuralFault,FaultHat)

%% Norm Calculation 

% Error = Error(:,end-3000:end);
% normError = [norm(Error(1,:))
%                     norm(Error(2,:))
%                     norm(Error(3,:))
%                     norm(Error(4,:))
%                     norm(Error(5,:))
%                     norm(Error(6,:))];
% 
% 
% uThrust = m*sqrt(u(1,:).^2+u(2,:).^2+(g+u(3,:)).^2);
% norm_uthrust = norm(uThrust); 
% norm_uphi = norm(u(4,:));
% norm_utheta = norm(u(5,:));
% norm_upsi = norm(u(6,:));
% 
% uNFTSMC = load('uNFTSMC_2');
% uNFTSMC = uNFTSMC.u;
% 
% xNFTSMC = load('xNFTSMC_2');
% xNFTSMC = xNFTSMC.x;
% 
% 
% uNFTSMC_1 = load('uNFTSMC');
% uNFTSMC_1 = uNFTSMC_1.u;
% 
% xNFTSMC_1 = load('xNFTSMC');
% xNFTSMC_1 = xNFTSMC_1.x;
% 
% 
% Error_NFTSMC = xNFTSMC(1:6,:)-XD;
% Error_NFTSMC = Error_NFTSMC(:,end-10000:end);
% normError_NFTSMC = [norm(Error_NFTSMC(1,:))
%                                     norm(Error_NFTSMC(2,:))
%                                     norm(Error_NFTSMC(3,:))
%                                     norm(Error_NFTSMC(4,:))
%                                     norm(Error_NFTSMC(5,:))
%                                     norm(Error_NFTSMC(6,:))];
% 
% uThrust_NFTSMC = m*sqrt(uNFTSMC(1,:).^2+uNFTSMC(2,:).^2+(g+uNFTSMC(3,:)).^2);
% norm_uthrust_NFTSMC = norm(uThrust_NFTSMC); 
% norm_uphi_NFTSMC = norm(uNFTSMC(4,:));
% norm_utheta_NFTSMC = norm(uNFTSMC(5,:));
% norm_upsi_NFTSMC = norm(uNFTSMC(6,:));
% 
% NormTable1 = [norm_uthrust           norm_uthrust_NFTSMC];
% 
% NormTable2 = [norm_uphi                norm_uphi_NFTSMC
%                        norm_utheta             norm_utheta_NFTSMC
%                        norm_upsi                norm_upsi_NFTSMC];


% Error_NN = FaultHat(:,end-10000:end)-SystemStructuralFault(:,end-10000:end);
% 

% normError_NN = [norm(Error_NN(1,:))
%                             norm(Error_NN(2,:))
%                             norm(Error_NN(3,:))
%                             norm(Error_NN(4,:))
%                             norm(Error_NN(5,:))
%                             norm(Error_NN(6,:))];
