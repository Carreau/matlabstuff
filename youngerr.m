function err = youngerr(params,Input,AO)
    %%    
    d0 = params(1);
    f0 = params(2);
    E = params(3);
    Drift = params(4);
    %%
    fitcurve = younghertz(Input,d0,f0,E,Drift);
    errvec = (fitcurve - AO) ;
    err = sum(errvec .^2) / length(AO) ;
    
    dmin = min(Input);
    dmax = max(Input);
    if( d0 < dmin )
      %err = err + 0.03*(dmin - d0)^2/err;
      %disp('malus d0 out of range');
    end
    if( d0 > dmax )
      %err = err + 0.03*(dmax - d0)^2/err;
      %disp('malus d0 out of range');
    end
    if (Drift >0 )
     %err = err+ 0.03*Drift^2/err ;
     %disp('malus positive drift');
    end
    if (E <0 )
     %disp('malus negative young');
     %err =  err - E*0.03/err
    end
    
    
end
