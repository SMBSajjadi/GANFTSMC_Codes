function J = GeneticFunction(K)

    Ts = 0.001;
    t = 0:Ts:50;
    T = pi/180;     % Conversion to Radians
    ALPHA = 20*T;
    GAMMA = 30*T;
    BETA = 10*T;
    FAULT_ANGLES = [ALPHA        BETA        GAMMA]';
    m = 1;
    g = 9.81;
    
    %% Initialization
    
    N = numel(t);
    n = 6;  % Number of DOF (X)
    
    x0 = [3 1 1 0 0 0 0 0 0 0 0 0];       %% Ics
    x = zeros(2*n,N);                          %% Compelete State Vector
    x(:,1) = x0;
    
    s = @(x) sin(x);
    c = @(x) cos(x);
    
    %% NFTSMC Parameters
    
    S = 0.1*ones(n,N);
    b = [0.01 0.01 0.01 0.1 0.1 0.1];       
    bPrime = [8 8 2 12 12 12];      
    Landa = 2;
    LandaPrime = 1.8;
    
    %% RBF Parameters

    GammaRBF = [0.01, 0.01,0.01,0.01,0.01,0.01]';
    nKernel = 10;           % Number of Kernel Functions
    bij = 0.1*ones(n,nKernel);
    cij = 0.01*ones(n,nKernel);
    
    %% Reference Signals
    
    CASE = 1;           %% Linear Trajectory
    % CASE = 2;       %% Eliptical Trajectory
    % CASE = 3;       %% BiLinear Trajectory
    [XD, XDotD, XDoubleDotD] = setDesiredTrajectory(t,CASE,n);
    
    phid = XD(4,:);
    thetad = XD(5,:);
    sayd = XD(6,:);
    
    %% Control Signal Initialization
    
    nU = 6;             %% iN cOnjunction with the Virtual Signals
    u = zeros(nU,N);
    u(3,1) = 0.1*m*g;
    Error = zeros(n,N);
    ErrorDot = 0.01*ones(size(Error));

    %% Fault Estimation Initialization
    
    FaultHat = 0.01*ones(n,N);
    wHatX = 0.01*ones(nKernel,N);
    wHatY = wHatX;
    wHatZ= wHatX;
    wHatPhi= wHatX;
    wHatTheta= wHatX;
    wHatPsi = wHatX;
    
    %% Main Loop
    
    for i=2:N  
        %% State Calculation
    
        [K1RK,wStar] =...
        Rotor2_DynamicForGenetic(t(i-1),...
        x(:,i-1),...
        u(:,i-1),...
        FAULT_ANGLES);
    
        x(:,i) = stateCalculation(K1RK,x(:,i-1),u(:,i-1),Ts,t(i-1),FAULT_ANGLES);
        
        %% Virtual Control Design
        
        thetad(i) = atan((u(1,i-1)*c(sayd(i))+u(2,i-1)*s(sayd(i)))/(g+u(3,i-1)));
        phid(i) = atan(c(thetad(i))*((u(1,i-1)*s(sayd(i))-u(2,i-1)*c(sayd(i)))/(g+u(3,i-1))));
        
        %% Desired Roll/Pittch Angles
        
        XD(4:5,i) = [phid(i), thetad(i)]'; 
    
        %% RBFNN-NFTSMC Design

        [u(:,i),S(:,i),Error(:,i),ErrorDot(:,i)] = ...
            RBFNNTSMC_GENETIC(wStar,x(:,i),n,XD(:,i),XDotD(:,i),XDoubleDotD(:,i),...
            b,bPrime, Landa,LandaPrime,K, FaultHat(:,i-1));
        
        %% RBF Estimator
        
        [wHatX(:,i),wHatY(:,i),wHatZ(:,i),wHatPhi(:,i),wHatTheta(:,i),wHatPsi(:,i), FaultHat(:,i)] = ...
        RBFNN_StructuralFaultEstimator(t,Error(:,i),ErrorDot(:,i),S(:,i),bij,cij,GammaRBF,...
        nKernel,wHatX(:,i-1),wHatY(:,i-1),wHatZ(:,i-1),...
        wHatPhi(:,i-1),wHatTheta(:,i-1),wHatPsi(:,i-1),bPrime,LandaPrime);

    end
   
    SUM_Error = sum(sum(Error.^2,2));
    Sum_u = sum(sum(u.^2,2));

    PenaltyU = 5;

    J = SUM_Error +...
            PenaltyU*Sum_u;

end