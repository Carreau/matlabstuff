function err = youngerr(params,Input,AO)
    %%    
    d0 = params(1);
    f0 = params(2);
    E = params(3);
    %%
    fitcurve = younghertz(Input,d0,f0,E);
    errvec = fitcurve - AO ;
    err = sum(errvec .^2);
end
