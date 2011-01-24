% fonction to reconvert raw datafile to more usables 
% resultsfiles, with more information 

function reconvert(patern)
[p,b,names] = dirr('.',patern,'name');
for i=1:length(names)
    nm = char(names(i));
    f = load(nm);
    newname = [nm(1:end-4) '_matthias.mat']
    save('newname','-struct','f', 'bead_radius', 'want_filter', 'f_3dB_filter', 'rate', 'T', 'l', 'f_lim_up', 'f_lim_low', 'nb', 'slopes', 'Traps', 'remove_peaks', 'eta_x', 'kappa_x', 'gamma_x', 'D_x', 'eta_y', 'kappa_y', 'gamma_y', 'D_y');
    clear f;
end

%'x', 'y', 'bead_radius', 'want_filter', 'f_3dB_filter', 'rate', 'T', 'l', 'f_lim_up', 'f_lim_low', 'nb', 'slopes', 'Traps', 'remove_peaks', 'eta_x', 'kappa_x', 'gamma_x', 'D_x', 'eta_y', 'kappa_y', 'gamma_y', 'D_y', 'Pbx', 'Pby', 'Pfitx', 'Pfity', 
