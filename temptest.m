function [s]=temptest( f )
    %tmp = load(str);
    %f  = tmp.save_data; 
    disp('********');
    disp('merci de donner le date et l''heure, ainsi qu''une description');
    disp('********');
    for i=1:length(f)
        s(i)=twoBeadApprochExperiment(f(i));
    end