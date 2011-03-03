function [young distance] = youngdistance(s)
    experiment = s.rawdata;
    start= experiment.touch_point;
    stop = experiment.appr_stop;
    distance = [experiment.d(start+49:stop)];
    young    = [experiment.E(50:end)];
    figure(1);
    hold on;
    plot(distance,young,'ro');
    xlabel('distance');
    ylabel('young modulus');
    title('young modulus vs distance');
    hold off;
    figure(2)
    xlabel('temps');
    ylabel('touch distance');
    hold on;
    plot(s.time_m,experiment.d(experiment.touch_point),'b+');
    hold off;
