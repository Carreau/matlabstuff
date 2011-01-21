function [eta_x_f,eta_y_f] = etafilter(eta_x,eta_y,r,min,max)
n = length(eta_x);
j=0;
for i=1:n
    ex = eta_x(i);
    ey = eta_y(i);
    rp = ex/ey;
    if( ex < max && ey < max && ex > min && ey > min && rp > r && rp < 1./r )
        j=j+1;
        eta_x_f(j)=ex;
        eta_y_f(j)=ey;
    end
end
