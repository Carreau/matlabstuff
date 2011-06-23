function err = ye(params,Input,AO)
    %%    
    d0 = params(1);
    E = params(2);
    %%
    fitcurve = yz(Input,d0,E);
    errvec = (fitcurve - AO) ;
    err = sum(errvec .^2) / length(AO) ;
    
end
