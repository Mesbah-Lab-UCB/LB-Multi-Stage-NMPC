function PreX = Pre(Omega,Ups,sys)

% Get sizes
nx = size(sys.A,2);
nu = size(sys.B,2);
nw = size(sys.E,2);

% Step 1: Calculate the projection of \Upsilon onto (x,u)
Y = [];
for i = 1:length(Ups)
    Y = [Y, projection(Ups(i),[1:nx+nu])];
end

% Step 2: Calculate the inverse map
Phi = [];
for i = 1:length(Ups)
    for j = 1:length(Omega)
        Aij = Ups(i).A;
        bij = Ups(i).b;
        Aij = [Aij ; Omega(j).A*[sys.A, sys.B, sys.E]];
        bij = [bij ; Omega(j).b];
        Phi = [Phi, Polyhedron('A',Aij,'b',bij)];
    end
end

% Step 3: Calculate the set difference
Ek = {};
for i = 1:length(Ups)
    Ek{i} = Ups(i);
end
for k = 1:length(Phi)
    for j = 1:length(Ek)
        if isEmptySet(Ek{j} & Phi(k))
            Ek{j} = Ek{j};
        else
            Ek{j} = Ek{j}\Phi(k);
        end
        Ek{j}.computeVRep;
    end
end

% Step 4: Compute the projection
Psi = [];
for i = 1:length(Ek)
    Psi = [Psi, projection(Ek{i},[1:nx+nu])];    
end

% Step 5: Compute the set difference
Sigma = Y\Psi;
Sigma.computeVRep;

% Step 6: Compute the projection
PreX = projection(Sigma,[1:nx]);
PreX.computeVRep;

end
