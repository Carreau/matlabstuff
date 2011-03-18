function [s]=temptest( f )
    %tmp = load(str);
    %f  = tmp.save_data; 
    disp('********');
    disp('merci de donner le date et l''heure, ainsi qu''une description');
    disp('********');
    for i=1:length(f)
        s(i)=twoBeadApprochExperiment(f(i));
    end
    disp('merci d''entrer la date et l''heure de mélange de l''actine');
    annee = input('annee : ');
    mois = input('mois : ');
    jour = input('jour : ');
    heures = input('heures : ');
    minutes = input('minutes : ');
    secondes = input('secondes : ');
    s.setDatevec(annee,mois,jour,heures,minutes,secondes);
    disp('veuillez entrer un concentration en Arp et cp');
    arp = input('arp : ');
    cp = input('cp : ');
    s.setArp(arp);
    s.setCp(cp);
end
