function [s]=temptest( str )
    tmp = load(str);
    f  = tmp.save_data;
    disp('cocuou');
    %disp('********');
    %disp('merci de donner le date et l''heure, ainsi qu''une description');
    %disp('********');
    for i=1:length(f)
        s(i)=twoBeadApprochExperiment(f(i));
    end
    matfilepath = [str(1:end-4),'_meta.mat'];
    var m;
    if(exist(matfilepath,'file'))
        disp('found meta data, load them');
        m = load(matfilepath);
    else
        disp('merci d''entrer la date et l''heure de mélange de l''actine');
        m.annee = input('annee : ');
        m.mois = input('mois : ');
        m.jour = input('jour : ');
        m.heures = input('heures : ');
        m.minutes = input('minutes : ');
        m.secondes = input('secondes : ');
        disp('veuillez entrer un concentration en Arp et cp');
        m.arp = input('arp : ');
        m.cp = input('cp : ');
        save(matfilepath,'-struct','m');
    end
    
    s.setDatevec(m.annee,m.mois,m.jour,m.heures,m.minutes,m.secondes);
    
    s.setArp(m.arp);
    s.setCp(m.cp);
end
