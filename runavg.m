function f = runavg(x,w)
    f = filter(ones(1,w)/w,1,x);
end
