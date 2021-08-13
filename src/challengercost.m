function J = challengercost(Gchall,y,b,iW)
%    function J = challengercost(Gchall,y,b,iW)

    ytilde = Gchall * b;
    ntilde = ytilde - y;
    J = ntilde'*iW*ntilde;
end

