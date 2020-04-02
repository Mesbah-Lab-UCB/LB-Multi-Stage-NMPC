%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Summary: Looks at a plant-model mismatch representation of disturbance
% w(x) = f(x,u) - a*x - b*u. In this case, w(x1) = 1-sin(x1.*1/(pi))
% Construct a GP regression from limited data
% Then build maximal robust control invariant set

% Written by:   Joel Paulson
% Date:         02/14/20
% Edited by:    Angelo D. Bonzanini
% Updated:      Feb 25 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize workspace

clearvars
rng(100,'twister')
uqlab


%% Specify problem data

% system matrices
A = [1, 1 ; 0, 1];
B = [1 ; 1];
E = [1 ; 1];
E = [1 ; 0];

% system sizes
nx = size(A,1);
nu = size(B,2);
nw = size(E,2);

% constraints
x_min = [-10 ; -10];
x_max = [ 10 ;  10];
u_min = -5;
u_max =  5;
X = Polyhedron('lb',x_min,'ub',x_max);
U = Polyhedron('lb',u_min,'ub',u_max);

% number of boxes
Nbox = 20;


%% Learn Guassian process (GP) model

% define model
ModelOpts.mString = '2-2*cos(X.*1.6/(pi))';
ModelOpts.isVectorized = true;
myModel = uq_createModel(ModelOpts);

% create input object
InputOpts.Marginals.Type = 'Uniform';
InputOpts.Marginals.Parameters = [x_min(1), x_max(1)];
myInput = uq_createInput(InputOpts);

% draw samples
Xsamp = uq_getSample(10,'MC');
Ysamp = uq_evalModel(myModel, Xsamp);

% train GP
MetaOpts.Type = 'Metamodel';
MetaOpts.MetaType = 'Kriging';
MetaOpts.ExpDesign.X = Xsamp;
MetaOpts.ExpDesign.Y = Ysamp;
myKrigingMat = uq_createModel(MetaOpts);


%% Construct bounding box

% loop over grid and find min and max values in each region
Xgrid = linspace(x_min(1),x_max(1),Nbox+1);
Delta = [];
for i = 1:Nbox
    % Define x values in selected window
    xU = Xgrid(i+1);
    xL = Xgrid(i);
    Xval = linspace(xL,xU,101)';
    % GP evaluation
    [YMeanMat,YVarMat] = uq_evalModel(myKrigingMat,Xval);
    points = [Xval, YMeanMat+sqrt(YVarMat)*norminv(1-0.01/2) ; Xval, YMeanMat-sqrt(YVarMat)*norminv(1-0.01/2)];
    % For each region find the convex hull and create a polyhedron object
    K = convhull(points);
    P = Polyhedron('V',points(K,:)); 
    Delta = [Delta, P.outerApprox()]; % Append outer approx (convert arbitrary shape to a state-dependent bounding box)
end

% find worst case disturbance values
Wv = projection(Delta,2);  % project on the w-axis (2nd dimension in this case)
Vlist = [];
for i = 1:Nbox
    Vlist = [Vlist ; Wv(i).V];
end

% Set in which all the uncertainties lie (worst-case bounds)
W = Polyhedron('lb',min(Vlist),'ub',max(Vlist));


%% Compute RCI set assuming bounding box

% Calculate the RCI set assuming outer-bounding box for W
sys = ULTISystem('A',A,'B',B,'E',E);
sys.x.min = x_min;
sys.x.max = x_max;
sys.u.min = u_min;
sys.u.max = u_max;
sys.d.min = min(W.V);
sys.d.max = max(W.V);

% RCI outer bound (with worst-case uncertainty)
Cinf_ob = sys.invariantSet('maxIterations',50);


%% Calculate RCI set using GP-based uncertainty description

% Define \Upsilon set when everything is independent
% = {(x,u,w) | x \in X, u \in U, w \in W}
Aup = blkdiag(X.A, U.A, W.A);
bup = [X.b ; U.b ; W.b];
Ups_convex = Polyhedron('A',Aup,'b',bup);

% Compute polytope for each convex cover set 
Ups = [];
for i = 1:length(Delta)
    % Need to add two columns of zeros becasue Upsilon refers to (x,u,w) = (x1, x2, u, w), whereas Delta refers to (x1,w)
    Ups_i = Polyhedron('A',[Aup ; [Delta(i).A(:,1), zeros(length(Delta(i).b),2), Delta(i).A(:,2)]],'b',[bup ; Delta(i).b]);
    Ups = [Ups, Ups_i];
end

% Backward reachability
N = 10;
Xlist{1} = Cinf_ob;
for k = 1:N
    % Print iteration
    fprintf('Iteration %g...',k)
    
    % Calculate the predecessor set
    Xnext = Pre(Xlist{k}, Ups, sys);

    % If set if convex, take the hull
    PU = PolyUnion(Xnext);
    if PU.isConvex()
        Xnext = PU.convexHull;
    end
    
    % Add to list
    Xlist{k+1} = Xnext;
    
    % Terminate if sets equal
    if Xlist{k+1} == Xlist{k}
        fprintf('converged\n')
        break
    end
    
    % Print complete statement
    fprintf('done\n')
end
Cinf = Xlist{k+1};


%% Select input that keeps you inside set

% pre-compute sets
Cinf_PU = PolyUnion(Cinf);
Cinf_CH = convexHull(Cinf_PU);
Delta_X1 = projection(Delta,1);

% big-M constant
M = 1000;

% filter constant
lambda_f = 0.5;

% output
C = [1, 1];

% target value
ytar = 5;

% run simulation
Nsim = 25;
Xsim = zeros(nx,Nsim+1);
Usim = zeros(nu,Nsim);
Wsim = zeros(nw,Nsim);
What = zeros(nw,Nsim);
Xsim(:,1) = Cinf(3).V(end,:)';

%{
for k = 1:Nsim
    % print start
    tic
    fprintf('simulation step %g / %g...', k, Nsim)
    
    % find the set W(x) for current state
    index = Delta_X1.contains(Xsim(1,k));
    index = find(index == 1);
    
    % Project onto the w-axis. If Xsim is contained in more than one boxes,
    % consider the box with the largest size (volume)
    if length(index) > 1
        Wx1 = projection(Delta(index(1)),2);
        Wx2 = projection(Delta(index(2)),2);
        if volume(Wx1) < volume(Wx2)
            Wx = Wx1;
        else
            Wx = Wx2;
        end
    else
        Wx = projection(Delta(index),2);
    end
    Wx.computeHRep;
    
    % calculate Pontryagin difference (erosion) for polygons (collection of polytopes)
    EWx = E*Wx;
    Cinf_next = (setMinus(X,(EWx)))\((X\Cinf) + (-EWx));
    npoly = length(Cinf_next);

%     % compute convex hull to improve bounding
%     Cinf_next_CH = convexHull(PolyUnion(Cinf_next));

    % define variables
    xss = sdpvar(nx,1);
    uss = sdpvar(nu,1);
    yss = sdpvar(nu,1);
    u = sdpvar(nu,1);
    d1 = binvar(npoly,1);
    
    % define optimization problem
    cons = [];
    cons = [cons, xss == A*xss + B*uss + E*What(:,k)];
    cons = [cons, C*xss == yss];
    cons = [cons, sum(d1) == 1];
    cons = [cons, U.A*u <= U.b];
    for j = 1:npoly
        cons = [cons, Cinf_next(j).A*(A*Xsim(:,k)+B*u) - Cinf_next(j).b <= M*(1-d1(j))];
    end
%     cons = [cons, Cinf_next_CH.A*(A*Xsim(:,k)+B*u) <= Cinf_next_CH.b];
    xnext = A*Xsim(:,k)+B*u+E*What(:,k);
    obj = (C*xnext-yss)'*(C*xnext-yss) + 100*(yss-ytar)^2;
    
    % solve optimization for current input
    ops = sdpsettings('solver','cplex','verbose',0,'debug',1);
    optimize(cons,obj,ops);
    
    % give values to plant
    Usim(:,k) = value(u);
    Wsim(:,k) = uq_evalModel(Xsim(1,k));
    Xsim(:,k+1) = A*Xsim(:,k) + B*Usim(:,k) + E*Wsim(:,k);
    
    % update estimate of disturbance
    What(:,k+1) = lambda_f*What(:,k) + (1-lambda_f)*Wsim(:,k);
    
    % print end message
    fprintf('took %g seconds\n', toc)
end


%% Plots
%}
% overlay maximal RCI calculated with w(x) and outerbounded w(x)
figure; hold on;
Cinf.plot('color','red','linestyle','-','edgecolor','red','linewidth',2)
Cinf_ob.plot('color','blue','linestyle','-','edgecolor','blue','linewidth',2)

% confidence regions and realized values of w(x1)
%%
figure; hold on;
Delta.plot('alpha',0.5)
Xval = linspace(sys.x.min(1),sys.x.max(1),1e3)';
Yval = uq_evalModel(myModel,Xval);
[YMeanMat,YVarMat] = uq_evalModel(myKrigingMat,Xval);
plot(Xval,YMeanMat, '--k', 'linewidth', 2);
plot(Xval, YMeanMat+sqrt(YVarMat)*norminv(1-0.01/2), '-r', 'linewidth', 2)
plot(Xval, YMeanMat-sqrt(YVarMat)*norminv(1-0.01/2), '-r', 'linewidth', 2)
xlabel('$x_1$', 'Interpreter', 'LaTeX')
ylabel('$\tilde{w}(x_1)$', 'Interpreter', 'LaTeX')
set(gca,'FontSize',15)
box on
% scatter(Xsim(1,1:end-1),Wsim(1,:),75,'bo','filled')
%%
%{
% phase plot
figure; hold on;
Cinf.plot('color','red','linestyle','-','edgecolor','red','linewidth',2)
plot(Xsim(1,:),Xsim(2,:),'-ko','linewidth',2);

% output plot
figure; hold on;
plot([0,Nsim],[ytar,ytar],'--r','linewidth',2);
plot(0:Nsim,C*Xsim,'-.ks','linewidth',2,'markersize',10);
%}
%%
%{
% Save mat file with relevant constraints
X_H = X.H;
U_H = U.H;
Cinf_CH_H = convexHull(PolyUnion(Cinf)).H;
Cinf_ob_H = Cinf_ob.H; % outer bound
Cinf_H = cell(size(Cinf,2),1);
for j=1:size(Cinf,2)
    Cinf_H{j,1} = Cinf(j).H;
end
%}

saveStr = ['constraintSetsM_', num2str(Nbox), '.mat'];
save(saveStr, 'Cinf', 'Cinf_CH', 'Cinf_ob', 'sys','X', 'U', 'myModel', 'myInput', 'myKrigingMat', 'Delta', 'Delta_X1')

