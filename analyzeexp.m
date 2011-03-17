function anaylseexp(s)
    fv = s.applyfun(@(x) x.fitvalue);
    
    % plot young modulus
    figure(1);
    plot([s.time_m],[fv.E],'+');
    title(s(1).commentaire);
    xlabel('time minutes');
    ylabel('young modulus');

    %plot touch distance
    figure(2);
    plot([s.time_m],[fv.d0],'+');
    xlabel('time minutes');
    ylabel('touch distance');
    title(s(1).commentaire);

    %plot error of fit
    figure(3);
    plot([s.time_m],[fv.err],'+');
    xlabel('time minutes');
    ylabel('fit error');
    title(s(1).commentaire);
end
