% function ment to only extract read files 
% and re-write a file containing only the force/distance
% and some information to link file together. 

function m=extract_force_distance(filename)
    %% load data
    %disp('chargement...');
    %g =    temptest('30_mars_10cp_25apr.mat');
    %g = [g temptest('30_mars_30cp_25apr.mat')];
    %g = [g temptest('30_mars_50cp_25apr.mat')];
    %disp('...');
    %g = [g temptest('8mars25Arp50cp.mat')];
    %g = [g temptest('8mars25arp00cp.mat')];
    disp('...');
    g= temptest(filename)
    %g = [g temptest('run1_1-mars-2011_25arp-10cp.mat')];
    %g = [g temptest('run3_1-mars-2011_25arp-30cp.mat')];
    %g = [g temptest('run4_1-mars-2011_25arp-30cp.mat')];

   
    m=[];
    for i=1:length(g)
       m(i).UUID     = java.util.UUID.randomUUID.toString;%#ok<AGROW>
       m(i).d        = g(i).bead_distance ;%#ok<AGROW>
       m(i).start    = g(i).moving_trap.event.appr.start;%#ok<AGROW>
       m(i).stop     = g(i).moving_trap.event.appr.stop;%#ok<AGROW>
       m(i).back     = g(i).moving_trap.event.retr_start;%#ok<AGROW>
       m(i).f        = g(i).still_trap_force.tangent;%#ok<AGROW>
       m(i).arp      = g(i).arp;%#ok<AGROW>
       m(i).cp       = g(i).cp;%#ok<AGROW>
       m(i).time_m   = g(i).time_m;%#ok<AGROW>
       m(i).parameters   = g(i).rawdata.parameters;%#ok<AGROW>
    end
    disp('done')


