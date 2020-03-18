%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary: Implements a learning-based multi-stage NMPC algorithm with
% state-dependent uncertainty. The plant-model mismatch is predicted online 
% using GP regression. Feasibility is guaranteed by enforcing an RCI
% constraint set, which is calculated in the main_gp_rci_v2.m script.

% Written by:   Angelo D. Bonzanini
% Date:         Feb 25 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize workspace

clear all
rng(100,'twister')
uqlab
import casadi.*

% Plot settings
Fontsize = 15;
Lwidth = 2;



%% Switch functionalities on/off
gpSwitch = 1;            % GP correction on/off
OFswitch = 0;            % Offset-free approach on/off
TrainOnly = 0;           % Carry out training only (when validating/testing)
saveSwitch = 1;          % Save outputs on/off
worstCase = 0;
useProj = 1;



%% Load outputs from main_gp_rci_v2.m script
% RCIout = load('constraintSetsM');    %15 boxes
% RCIout = load('constraintSetsM_20'); %20 boxes
RCIout = load('constraintSetsM_30'); %30 boxes

% Extract relevant objects
sys = RCIout.sys;
X = RCIout.X;
U = RCIout.U;
Cinf = RCIout.Cinf;
Cinf_CH = RCIout.Cinf_CH;
Cinf_ob = RCIout.Cinf_ob;
myModel = RCIout.myModel;
myInput = RCIout.myInput;
Delta = RCIout.Delta;
Delta_X1 = projection(Delta,1);


nx=size(sys.A,2);
nu=size(sys.B,2);
% Double integrator model
A = sys.A;
B = sys.B;
C = eye(2,2);
E = sys.E;

%% Learn Guassian process (GP) model

% draw samples
Xsamp = uq_getSample(myInput,8,'MC');
Ysamp = uq_evalModel(myModel, Xsamp);

% train GP
MetaOpts.Type = 'Metamodel';
MetaOpts.Scaling=0;
MetaOpts.MetaType = 'Kriging';
% MetaOpts.Corr.Family = 'matern-3_2';
MetaOpts.Corr.Family = 'Gaussian';
MetaOpts.Corr.Isotropic = 1;

MetaOpts.ExpDesign.X = Xsamp;
MetaOpts.ExpDesign.Y = Ysamp;
myKrigingMat = uq_createModel(MetaOpts);

% Extract hyperparameters
theta = myKrigingMat.Internal.Kriging.Optim.Theta;
sigmaSQ = myKrigingMat.Kriging.sigmaSQ;
b = myKrigingMat.Kriging.beta;
Rmat = myKrigingMat.Internal.Kriging.GP.R;
Fmat = myKrigingMat.Internal.Kriging.Trend.F;

% Test Data
Xtest = linspace(-10, 10, 20)';
[Ypred,Yvar] = uq_evalModel(myKrigingMat,Xtest);
Ytest = uq_evalModel(myModel, Xtest);

% Calculate kernel "manually" for consistency check
XXsamp = myKrigingMat.ExpDesign.U;
Rmanual = kernelFn(Xsamp, Xsamp, theta);
if abs(Rmat-Rmanual)<=1e-4
    fprintf('\nKernel consistency test passed!\n')
else
    error('Kernel consistency test failed!')
end

% Plot test data and predictions
figure(1)
hold on
h1 = plot(Xtest, Ytest, 'b.', 'MarkerSize', 10);
h2 = plot(Xtest, Ypred, 'r-', 'MarkerSize', 10);
h3 = plot(Xtest, Ypred+3*sqrt(Yvar), 'k--');
plot(Xtest, Ypred-3*sqrt(Yvar), 'k--')
xlabel('X_{test}')
ylabel('Y')
set(gca,'FontSize',Fontsize)
box on

sigma = max(abs(3*sqrt(Yvar)));

if TrainOnly==1
    return
end



%% Project into maximal robust control invariant set
%{
if useProj==1
     [explicit_controller, mptsol,diagn,Z,Valuefcn,Optimizer] = CinfProjection_ob(X, U, Cinf, Cinf_ob, Delta, Delta_X1, sys, worstCase);
     save('uexp.mat', 'explicit_controller')
%     load('uexp')
end
%}






%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MULTI-STAGE NMPC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MPC Parameters

% Define cost parameters
Q = [1, 0; 0, 1];
R = 0.1;
PN = Q;                            % Terminal cost weight

nx = size(Q,1);
nu = size(R,1);

Np = 4;                             % Prediction horizon
N = 20;                             % Simulation horizon
N_robust = 2;                       % Robust horizon for multistage MPC

% Initial point
yi = [-10;-4.5];

% Set point(s) and time(s) at which the reference changes
ysp1 = [10;0];
tChange = 10;
ysp2 = [0;0];
tChange2 = 9999;
ysp3 = [0;0];


yss = ysp1;
uss = 0;





%% Variables and model
% Declare the variables
x1 = MX.sym('x1');
x2 = MX.sym('x2');
u1 = MX.sym('u1');
x = [x1;x2];
u = u1;

ss = MX.sym('xssVar', 3);
wNoise = MX.sym('wNoise', 2);


% Double integrator model
A = sys.A;
B = sys.B;
C = eye(2,2);


%% Offset-free and estimation matrices (if needed)

%Offset-free tracking parameters
Bd = eye(2,2);
Cd = zeros(2,2);
Haug = eye(2,2);
Aaug = [eye(2,2)-A, -B; Haug*C, zeros(nx, nu)];


% Steady-state Kalman Gain
Adist = [A, Bd; zeros(2,2), eye(2,2)];
Bdist = [B;zeros(nx, nu)];
Cdist = [C, Cd];
Qnoise = 10*diag([0.1, 0.1, 0.1, 0.1]);
Rnoise = 10*diag([0.01, 0.01]);
[Pinf, ~, Kinf] = dare(Adist', Cdist', Qnoise, Rnoise);
Kinf = -Kinf;
% LQR gain
[Plqr, ~, Klqr] = dare(A, B, Q, R);
Klqr = -Klqr;


%% Solve the setpoint calculator multiparametrically to reduce computation time
mpSP = spCalculator(Aaug, Bd, Cd, X, U, Haug, sys);



%% Constraints
% Already loaded as X and U
% RCI set already loaded as Cinf
% RCI outer bound set already loaded as Cinf_ob




%% Define dynamics and cost as Casadi functions

% Dynamics and stage cost
xNext = A*x+B*u + wNoise;
y = C*x;
Lstage = (y-ss(1:2))'*Q*(y-ss(1:2)) + (u-ss(3:end))'*R*(u-ss(3:end));

% Functions for nominal vs. real dynamics
F = Function('F', {x,u,wNoise,ss}, {xNext,Lstage},{'x','u', 'wNoise', 'ss'},{'xNext', 'Lstage'});


% Variable inputs for test data (since these will change within the loop)
xTest = MX.sym('xTest', 1,1);


% CadADi function for GP prediction
Ks = kernelFn(Xsamp, xTest, theta);
Kss = kernelFn(xTest, xTest, theta);
Kinv = inv(Rmat);
Var = Kss - Ks'*Kinv*Ks;



if worstCase==0 && gpSwitch==1
    Yout = b+Ks'*Kinv*(Ysamp-b);
elseif worstCase==1 && gpSwitch==1
    Yout = max(abs(Ysamp));
elseif gpSwitch==0
    Yout = zeros(1,1);
else
    error('INVALID CASE STUDY COMBINATION!!!')
end

Fpred = Function('Fpred', {xTest}, {Yout},{'xTest'},{'Yout'});
Fvar = Function('Fvar', {xTest}, {Var},{'xTest'},{'Yout'});


% Confirm that the manually coded kernel yields the same results
Ytest2 = Fpred(Xtest);
figure(1)
h4 = plot(Xtest, full(Ytest2), 'mx');
legend([h1, h2, h3, h4], 'Test Data', 'GP Predictions', '99% Confidence Interval', 'Manual Kernel')



%% Initialize MPC

% Scenarios
scenario_idx = [1, 0, -1];     % Multiplier of the additive w(x,u) of the GP in each scenario

% Build the scenario matrix with all the combinations. Dimension (Ncases^Nrobust, Nrobust)
scenario_mat = combvec(scenario_idx, scenario_idx)';
N_scenarios = length(scenario_mat);

% Weights for cost function
w_i = (1/N_scenarios)*ones(N_scenarios,1);


% Initialize vectors to store the predicted trajectories for each scenario
y1S = zeros(N_scenarios, Np+1);
y2S = zeros(N_scenarios, Np+1);
u1S = zeros(N_scenarios, Np);
u2S = zeros(N_scenarios, Np);
gp1S = zeros(N_scenarios, Np);
gp2S = zeros(N_scenarios, Np);

% Initialize vectors to store the predicted inputs and outputs for one OCP loop
uopt = zeros(nu, Np);
yModel = zeros(nx, N+1);

% Initialize vectors to store the actual applied inputs and measured outputs
uOptSeq = zeros(nu, N);
fopt = zeros(N,1);
yTr = zeros(nx, N+1);


% Initialization
ssPlot = [ysp1(1);ysp1(2)];
ssPlot = [ssPlot, ssPlot];
YcMat = [];
dhat = zeros(2,1);
xki = yi;
xhati = xki;
yTr(:, 1)= yi;



%% MPC Loop
Tstart = tic;

for k = 1:N
    xki = xhati;
    Jactual = 0;
    gp1_opt = zeros(1, Np+1);
    gp2_opt = zeros(1, Np+1);
    
    
%     fprintf('\n\n################################# NEW OPTIMIZATION #################################\n\n')
    
    %   At each step k the entire OCP for all scenarios is solved!
    %   Therefore, we need to intialize empty variables for each
    %   step k.
    
    % Start with an empty NLP
    w=[];    %Array of all the variables we will be optimizing over
    w0 = [];
    lbw = [];
    ubw = [];
    discrete = [];
    J = 0;
    g=[];
    lbg = [];
    ubg = [];
    strMat = {};
    
    % Check for errors
    checkSc = [];
    feasMeasure = cell(N,1);
%     disp(xhati) 
    
    
    % MPC LOOP FOR DIFFERENT SCENARIOS - ROBUST HORIZON = 2
    for n_sc =1:length(scenario_mat)
        sc_vec = scenario_mat(n_sc, :)';
        
        rng(n_sc)
       
        wReal = [0*normrnd(0., 0.3, [1,N+1]);0*normrnd(0., 0.3, [1,N+1])];

        %     "Lift" initial conditions. Note that the initial node
        %     is the same for all scenarios, so the double index is not
        %     required.
        
        Xk = MX.sym(char(join(["X0","_",string(n_sc)])), nx);
        strMat = [strMat;{char(join(["X0","_",string(n_sc)]))}];
        strMat = [strMat;strMat{end}];
        w = [w;Xk];
        lbw = [lbw;xhati];
        ubw = [ubw;xhati];
        w0 = [w0;zeros(nx,1)];
        discrete =[discrete;zeros(nx,1)];
        
        Yk = MX.sym(char(join(['Y0','_',string(n_sc)])), nx);
        strMat = [strMat;{char(join(['Y0','_',string(n_sc)]))}];
        strMat = [strMat;strMat{end}];
        w = [w;Yk];
        lbw = [lbw; -inf*ones(nx,1)];
        ubw = [ubw; inf*ones(nx,1)];
        w0 = [w0; zeros(nx,1)];
        discrete =[discrete; zeros(nx,1)];
        
        
        sdGP = [0;0]; %(max(abs(sigma))*[1;0]
        wPred = [Fpred(xki(1)), 0];
%         wPred = [wFn(xki(1)), 0];
        
        
        % Optimal Control Problem - Open Loop Optimization
        
        % Formulate the NLP
        for i = 1:Np

            % New NLP variable for the control
            Uk = MX.sym(char(join(['U_',string(i),'_',string(n_sc)])), nu);
            strMat = [strMat;{char(join(['U_',string(i),'_',string(n_sc)]))}];
            w   = [w;Uk];
            lbw = [lbw; U.V(2)];
            ubw = [ubw;U.V(1)];
            w0  = [w0; zeros(nu,1)];
            discrete =[discrete;zeros(nu,1)];
            
            
            YGP = MX.sym(char(join(['YGP_',string(i),'_',string(n_sc)])), nx);
            strMat = [strMat;{char(join(['YGP_',string(i),'_',string(n_sc)]))}];
            strMat = [strMat;strMat{end}];
            w   = [w;YGP];
            lbw = [lbw; -9999*ones(nx,1)];
            ubw = [ubw; 9999*ones(nx,1)];
            w0  = [w0;zeros(nx,1)];
            discrete =[discrete;zeros(nx,1)];
            
            
            
            % Integrate until the end of the interval
            if i<=N_robust
                [Xk_end, Jstage] = F(Xk, Uk, gpSwitch*(YGP+sc_vec(i)*sdGP),[yss;uss]);
            else
                if worstCase==0
                    [Xk_end, Jstage] = F(Xk, Uk, gpSwitch*YGP,[yss;uss]);
                else
                    [Xk_end, Jstage] = F(Xk, Uk, [0;0],[yss;uss]);
                end
            end
                

            % Yk_end = mtimes(C, Xk_end)+0*YGP
            J=J+w_i(n_sc)*Jstage;
            % Penalize abrupt changes
            %J = J + mtimes(mtimes((Uk-uopt[:,i]).T, RR), Uk-uopt[:,i]) #+ mtimes(mtimes((Yk_end-Yk).T, QQ), Yk_end-Yk)
            
            
            % State-dependent uncertainty
            xx = Xk;
            if i<=Np
                wPred = [wPred;Fpred(Xk(1)), 0]; % Correct for first prediction   
            else
                wPred = [wPred; 0, 0]; % Correct for future predictions
            end
                
                

            g   = [g;Yk-C*Xk];
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
            
            g = [g;YGP-gpSwitch*wPred(i,:)'];
            lbg = [lbg;zeros(nx,1)];
            ubg = [ubg;zeros(nx,1)];
 
            % New NLP variable for state at end of interval
            Xk = MX.sym(char(join(['X_',string(i+1),'_',string(n_sc)])), nx);
            strMat = [strMat;{char(join(['X_',string(i+1),'_',string(n_sc)]))}];
            strMat = [strMat;strMat{end}];
            w   = [w;Xk];
            lbw = [lbw;-inf*ones(nx,1)];
            ubw = [ubw; inf*ones(nx,1)];
            w0  = [w0;zeros(nx,1)];
            discrete =[discrete;zeros(nx,1)];
            
            
            
            Yk = MX.sym(char(join(['Y_',string(i+1),'_', string(n_sc)])), nx);
            strMat = [strMat;{char(join(['Y_',string(i+1),'_', string(n_sc)]))}];
            strMat = [strMat;strMat{end}];
            w   = [w;Yk];
            ubw = [ubw;X.V(2,1);X.V(2,2)];
            lbw = [lbw; X.V(4,1);X.V(4,2)];
            w0  = [w0;zeros(nx,1)];
            discrete =[discrete;zeros(nx,1)];
            
            g   = [g;Xk_end-Xk];
            lbg = [lbg; zeros(nx,1)];
            ubg = [ubg; zeros(nx,1)];
            
            %{
            if i==2
            Ybin = MX.sym(char(join(['Ybin_',string(i),'_', string(n_sc)])), npoly,1);
            strMat = [strMat;{char(join(['Y_',string(i),'_', string(n_sc)]))}]; 
            w   = [w;Ybin];
            lbw = [lbw;zeros(npoly,1)];
            ubw = [ubw; ones(npoly,1)];
            w0  = [w0;zeros(npoly,1)];
            discrete =[discrete;ones(npoly,1)];
            sumCon = 0;
            M=50;
                for ii=1:npoly
                    strMat = [strMat;strMat{end}]; 
                    g   = [g;Cinf_next(ii).A*Xk-Cinf_next(ii).b-M*(1-Ybin(ii))];
                    lbg = [lbg; -inf*ones(length(Cinf_next(ii).b),1)];
                    ubg = [ubg; zeros(length(Cinf_next(ii).b),1)];
                    sumCon = sumCon + Ybin(ii);
                end
                g = [g; sumCon-1];
                lbg = [lbg;0];
                ubg = [ubg;0];
                
            end
            %}
            
            
        end
        % Terminal cost and constraints (Xk --> i+1)
        % Terminal Cost
        J = J + w_i(n_sc)*(Yk-yss)'*PN*(Yk-yss);
        
        % Equality constraint to make sure that Yk at the last step is equal to  C*Xk
        g = [g;Yk-C*Xk];
        lbg = [lbg;zeros(nx,1)];
        ubg = [ubg;zeros(nx,1)];    
        
        
    end
    



    
    conCheck0=[];
    conCheck1=[];    
    
    % First split - Find all the U_0's
    %%
    u0idx = [];
    u1idx = [];
    for l=1:length(strMat)
        st0= strfind(strMat{l}, 'U_ 1');
        st1= strfind(strMat{l}, 'U_ 2');
        if isempty(st0)==0
            u0idx =[u0idx;l];
        end
        if isempty(st1)==0
            u1idx =[u1idx;l];
        end
    end
    %%
    
    for con_idx = 1:length(u0idx(1:end-1))
        g = [g;w(u0idx(con_idx))-w(u0idx(con_idx+1))];
        conCheck0 = [conCheck0;w(u0idx(con_idx))-w(u0idx(con_idx+1))];
        lbg = [lbg;zeros(nu,1)];
        ubg = [ubg;zeros(nu,1)];
    end

    Nbranch = length(scenario_idx); % branches per node
    Nrepeat = length(w)/N_scenarios;
    %%
    
    % Second split
    
    % Group second scenarios based on the first split and then examine each
    % group of branches individually
    splitIdx = cell(1);
    for l = 1:length(scenario_idx)
        splitIdx{l,1} = find(scenario_mat(:,1)==scenario_idx(l));
    end
    
    % The scenarios U_2_splitIdx{l} should be equal, since they stem from the
    % same parent node
    
    for l = 1:length(scenario_idx)
        for ll =1:length(splitIdx{l})-1
            g = [g;w(u1idx(splitIdx{l}(ll)))-w(u1idx(splitIdx{l}(ll+1)))];
            conCheck1 = [conCheck1;w(u1idx(splitIdx{l}(ll)))-w(u1idx(splitIdx{l}(ll+1)))];
            lbg = [lbg;zeros(nu,1)];
            ubg = [ubg;zeros(nu,1)];
        end
    end
    
    
  
    %%
    % Create an NLP solver
    prob = struct('f', J, 'x', w, 'g',g);
%     sol_opts = struct('ipopt.print_level',0, 'ipopt.max_cpu_time',10, "discrete", discrete);
    sol_opts = struct('discrete', discrete);
    sol_opts.ipopt.max_iter = 500;

    solver = nlpsol('solver', 'ipopt', prob, sol_opts);
    
    % Solve the NLP
    sol = solver('x0', w0, 'lbx',lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

    w_opt = full(sol.x);
    J_opt = full(sol.f);
    
    
    % These only include the first scenario --> first case in scenario mat
    y1_opt = w_opt(3:7:7*Np+1)';
    y2_opt = w_opt(4:7:7*Np+1)';
    u1_opt = w_opt(5:7:7*Np+1)';
    gp1_opt = w_opt(6:7:7*Np+1)';
    gp2_opt = w_opt(7:7:7*Np+1)';
    
    feasMeasure{k} = solver.stats.return_status;
    
    
    % Perform plant simulations for the next step using the plant model
    uopt = u1_opt;
    
    % specify to use the projection or just the DNN
    if useProj == 1
        fprintf('simulation step %g / %g...', k, N)
        if worstCase==0
    %       uopt(:,1) = explicit_controller([xhati;uopt(:,1)]);
            uopt(:,1) = CinfProjection(xhati, uopt(:,1), X, U, Cinf, Delta, Delta_X1, sys);
        else
           Delta_X1_ob = PolyUnion(Delta_X1).outerApprox();
           Delta_ob = PolyUnion(Delta).outerApprox();
           uopt(:,1) = CinfProjection(xhati, uopt(:,1), X, U, Cinf_ob, Delta_ob, Delta_X1_ob, sys); 
        end
        fprintf('took %g seconds\n', toc)
    else
        uopt(:,1) = uopt(:,1);
    end
    
    uOptSeq(:,k) = uopt(:,1);
    fopt(k) = J_opt;
    
    %     Update intial condition. Can change the if statement to define 
    %     when the plant-model mismatch is introduced (e.g. glass-to-metal
    %     transition).
    
    dist = [wFn(xki(1)); 0];
    dhat = [Fpred(xki(1));0];
    
    dhat = dist;
    xki = A*xki+B*uopt(:,1) + dist +0*wReal(:,k);
    yki = C*xki;
    
    yTr(1,k+1) =yki(1);
    yTr(2,k+1) =yki(2);
    
    
    % State Feedback (comment out for output feedback)!!
    xhati = yki;
   

    %
    if k>=0 && k<tChange-1
        yspCase = ysp1;
    elseif k>=tChange-1 && k<tChange2-1
        yspCase = ysp2;
    elseif k>=tChange2-1
        yspCase = ysp3;
    end
        
    % Setpoint calculator 
%     sp_opt = spCalculator(yspCase, dhat, Aaug, Bd, Cd, X, U, Haug, sys);
    sp_opt = mpSP([yspCase;dhat]);
    yss = sp_opt(1:2);
    uss = sp_opt(3:end);
    ssPlot=[ssPlot, yspCase];    
    
end
%%

% Tend = toc;
% disp(['Total time = ', num2str(Tend)]);
% disp(['Average time = ', num2str(Tend/N)]);


figure(2)
subplot(2,1,1)
hold on
plot([0:N], yTr(1,:), 'LineWidth', Lwidth)
plot([0:N], ssPlot(1,1:end-1), 'k--', 'LineWidth', Lwidth)
ylim([-11, 11])
ylabel('x_1')
set(gca,'FontSize',Fontsize)
box on
subplot(2,1,2)
hold on
plot([0:N], yTr(2,:), 'LineWidth', Lwidth)
plot([0:N], ssPlot(2,1:end-1), 'k--', 'LineWidth', Lwidth)
ylim([-11, 11])
xlabel('Time Step (k)')
ylabel('x_2')
set(gca,'FontSize',Fontsize)
box on


figure(3)
hold on
plot(X, 'color', [1, 1, 1]*0.9, 'LineWidth', Lwidth, 'LineStyle', '--')
plot(yTr(1,:), yTr(2,:), 'LineWidth', Lwidth)
plot(ysp1(1), ysp1(2), 'kx', 'Markersize', 10)
plot(ysp2(1), ysp2(2), 'kx', 'Markersize', 10)

if saveSwitch==1
    save(['../Output-Data-Files/LB-MS-MPC_', datestr(now,'YYYY-mm-dd_HH_MM_SS'), ], 'X', 'Cinf', 'Cinf_ob', 'yTr')
end








