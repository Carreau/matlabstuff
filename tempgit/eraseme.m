x=[0:0.01:1];
f = @(alpha,x) (1/alpha).^(3/2).*((1-alpha.*x)./(1-x)).^(3/2)

j=1;
for i=[0.90:0.025:1.1]
   ff(j,:)=f(i,x); 
   j=j+1;
end


plot(x,ff)
axis([0    1.0000   -0.0000   3.0000]);
 
 
 