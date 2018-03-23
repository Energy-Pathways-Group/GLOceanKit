s = 212;
sD = s * exp(-5000/1300);
nuAxis = linspace(sD,s/2,100);

g = zeros(size(nuAxis));
for i=1:length(g)
    nu = nuAxis(i);
   g(i) = (((s/2).^(2*nu))/(2*(nu^3)*gamma(nu)^2)) * genHyper([nu, nu+1/2],[nu+1,nu+1,2*nu+1],-s*s);
end

figure, plot(nuAxis,g)