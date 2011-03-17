function err = youngerr(params,Input,AO)
    %%    
    d0 = params(1);
    f0 = params(2);
    E = params(3);
    Drift = params(4);
    %%
    fitcurve = younghertz(Input,d0,f0,E,Drift);
    errvec = fitcurve - AO ;
    err = sum(errvec .^2);
end
