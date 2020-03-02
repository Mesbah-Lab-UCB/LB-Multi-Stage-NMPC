function [explicit_controller, mptsol,diagn,Z,Valuefcn,Optimizer] = CinfProjection_ob(X, U, Cinf, Cinf_ob, Delta, Delta_X1, sys, worstCase)

A = sys.A;
B = sys.B;
E = sys.E;
nx = size(A,2);
nu = size(B, 2);

% Define problem to project into Cinf
xcurrent = sdpvar(nx,1);
uexplicit = sdpvar(nu,1);
uproject = sdpvar(nu,1);

constraints = [];
objective = 0;
xnext = A*xcurrent + B*uproject;



%     index = Delta_X1.contains(xnext);
%     index = find(index == 1);
%     Wx = projection(Delta(index),2);
Wx = PolyUnion(projection(Delta, 2)).outerApprox();
 %{
 M = 50;

 for l=1:size(Delta,1)
     constraints = [constraints, Delta(l).A*xnext<=Delta(l).b+M*(1-ybin(l))];
 end
 constraints = [constraints; sum(ybin)==1];    
%}

% calculate Pontryagin difference (erosion) for polygons (collection of polytopes)
EWx = E*Wx;
% Cinf_next = (setMinus(X,(EWx)))\((X\Cinf) + (-EWx));
% Wx = PolyUnion(Delta).outerApprox();
if worstCase==0
    Cinf_next = PolyUnion(Cinf).outerApprox() - EWx;
else
    Cinf_next = Cinf_ob-EWx
end
Cinf_next.computeVRep();
% Cinf_next = PolyUnion(Cinf_next).outerApprox();


constraints = [constraints, ismember(xnext, Cinf_next)];
constraints = [constraints, U.A*uproject <= U.b];
objective = objective + (uexplicit - uproject)'*(uexplicit - uproject);

% Add constraints on the explicit variable to bound the size of the mp map
% constraints = [constraints, -100*(sys.u.max-sys.u.min)' + sys.u.min' <= uexplicit <= sys.u.max' + 100*(sys.u.max-sys.u.min)'];

% Create optimizer object
ops = sdpsettings('verbose',0);
explicit_controller = optimizer(constraints,objective,ops,[xcurrent;uexplicit],[uproject]);

% Calculate the explicit solution using yalmip
[mptsol,diagn,Z,Valuefcn,Optimizer] = solvemp(constraints,objective ,ops,[xcurrent;uexplicit],[uproject]);
end

