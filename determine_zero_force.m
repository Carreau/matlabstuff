function [zero] = determine_zero_force(data)

s = data.start;
l= 10000;
dirforce = data.f(s:s+l)*1e12;
invforce = data.f(end-l:end)*1e12;
zero = mean((dirforce+invforce)/2)*1e-12;

end