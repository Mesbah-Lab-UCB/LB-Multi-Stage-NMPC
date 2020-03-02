function mpSP = spCalculator(Aaug, Bd, Cd, Ycon, Ucon, Haug, sys)

% Model
A = sys.A;
B = sys.B;
C = eye(2,2);

% Dimensions
nx = size(A,2);
nu = size(B,2);

% Objective function weights
Qs = [1, 0; 0, 0.1];

% Define variables
ysp = sdpvar(nx,1);
usp = sdpvar(nu,1);

yspNominal = sdpvar(nx,1);
dhat = sdpvar(nx,1);

% Objective function
objective = (ysp-yspNominal)'*Qs*(ysp-yspNominal);

constraints = [ysp==A*ysp+B*usp+dhat;
                Ycon.A*ysp<=Ycon.b;
                Ucon.A*usp<=Ucon.b];


% Create optimizer object
ops = sdpsettings('verbose',0);
mpSP = optimizer(constraints,objective,ops,[yspNominal;dhat],[ysp;usp]);

% Calculate the explicit solution using yalmip
solvemp(constraints,objective ,ops,[yspNominal;dhat],[ysp;usp]);


end

