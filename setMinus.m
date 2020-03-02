function Pdiff = setMinus(P, S)

% Subtract the support of S for each inequality of P
Hn = [P.A P.b-S.support(P.A')];
if any(Hn(:,end)==-Inf)
    % empty polyhedron in the same dimension
    Pdiff = Polyhedron.emptySet(P.Dim);
    return;
end
% remove remaining Inf rows
Hn((Hn(:,end))==Inf,:) = [];

Pdiff = Polyhedron('H', [Hn], 'He', P.He);

end
