function ret = yz(x,d0,E)
    %% calcul du module d'young 
    % Usage : 
    % younghertz(x,d0,f0,E)
    % F=4/3*E/(1-(nu)^2).*d.^(3/2)*sqrt(r)+f0;
    % avec 
    % F force sur la bille
    % E module d'young
    % nu coefficient de poisson 
    % d ?
    % ?
    % f0
    nu= 0.5 ; %
    r = 4.5;
    f=@(d) 4/3*E/(1-(nu)^2) .* d .^(3/2) *sqrt(r); 
    o = ones(size(x));
    f0=0;
    %pos = (x>d0*o)*(f0*o+f((x-d0*o))) ;
    ret = (f0*o+f(d0-x)) ;
    %neg = (x>d0) .* (f0);
    %ret = pos+neg;
end
