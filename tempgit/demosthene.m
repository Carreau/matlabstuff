function [temps,cell]=demosthene()
    temps = [];
    while(isempty(temps))
        temps = input('Num�ro de la colonne temps : ');
    end
    n=[];
    while(isempty(n))
        n=input('Num�ro de la colonne pour premi�re cellule : ');
    end
    i=2;
    cell = n;
    while(~isempty(n))
        n=input(sprintf('numero de la colonne pour la cellule %d, ou rien pour continuer : ',i+1));
        if(~isempty(n))
            cell(i)=n;
        end
        i=i+1;
    end
end