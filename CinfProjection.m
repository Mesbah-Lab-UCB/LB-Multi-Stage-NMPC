function uProj = CinfProjection(x, u, X, U, Cinf, Delta, Delta_X1, sys)

fprintf("\n Projecting to safe set...")
% Project optimal input onto maximal robust control invariant set

A = sys.A;
B = sys.B;
E = sys.E;
nx = size(sys.A,2);
nu = size(sys.B, 2);


% Define problem to project into Cinf
xcurrent = sdpvar(nx,1);
uproject = sdpvar(nu,1);
xnext = A*xcurrent + B*uproject;

% Initialize objective fn and constraints
constraints = [];
objective = 0;


% find the set W(x) for current state
index = Delta_X1.contains(min(10, x(1)));
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


%ismember constraint
constraints = [constraints, ismember(xnext, Cinf_next)];
% constraints = [constraints, ismember(xcurrent, Cinf)];
%{
for l=1:length(npoly)
    constraints = [constraints, Cinf_next(l).A*xnext <= Cinf_next(l).b+M*(1-ybinNext(l))];
    constraints = [constraints, Cinf(l).A*xcurrent <= Cinf(l).b+M*(1-ybin(l))];
end
%}
constraints = [constraints, U.A*uproject <= U.b];
objective = objective + (u - uproject)'*(u - uproject);

% Add constraints on the explicit variable to bound the size of the mp map
% constraints = [constraints, -100*(sys.u.max-sys.u.min)' + sys.u.min' <= uexplicit <= sys.u.max' + 100*(sys.u.max-sys.u.min)'];


% Create optimizer object
ops = sdpsettings('verbose',1);
optimize(constraints,objective,ops);
uProj = value(uproject);
%%


fprintf("done \n")
end

