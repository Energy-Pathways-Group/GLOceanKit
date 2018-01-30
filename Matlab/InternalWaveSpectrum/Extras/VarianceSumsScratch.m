j = 1:300;

j_star = 3;
H1 = (1+j).^(-5/2);
H1_norm = 1/sum(H1);
H1 = @(j) H1_norm*(1+j).^(-5/2);

sum(H1(j))

sum(H1(j/j_star)/j_star)

% return

H2 = (j+1).^(-5/2);
H2_norm = 1/sum(H2);
H2 = H2_norm*(j+1).^(-5/2);

sum(H2)

D = 5000;
m = j*pi/5000;
m_star = j_star*pi/5000;
Hm = (H1_norm/m_star)*(1+m/m_star).^(-5/2);

sum(Hm(3:end))*m(1)

% return

% figure
% plot(j,H1), hold on, plot(j,H2)

s=1; t=2.5;
A = s*gamma(t/s)/(gamma(1/s)*gamma((t-1)/s));
H = @(x) A*(1+x.^s).^(-t/s);

x = linspace(0,100,1000)';
trapz(x,H(x))
sum(H(j/j_star)/j_star)

sum(H(x))*(x(2)-x(1))

sum(H(j))

% Why doesn't this work?
% Because the discrete sum is too finite.
trapz(m,H(m/m_star)/m_star)
A = sum(H(m/m_star))*(m(2)-m(1))/m_star

figure
plot(x,H(x)), xlog, ylog, hold on
scatter(m/m_star,H(m/m_star)*(m(2)-m(1))/m_star/A)


j_star = 3;
H1 = (j_star+j).^(-5/2);
H1_norm = 1/sum(H1);
H1 = H1_norm*(1+j).^(-5/2);

dz = 1; % assume 1 meter scale
r = j_star;
p = 1.250;
A2 = 2*r^(2*p-1)*gamma(p)/gamma(p-0.5)/sqrt(pi);
c_alpha = gamma(0.5)*gamma(p-0.5)/gamma(p)/(2); % from Jonathan...this is correct.
A2 = r^(2*p-1)/c_alpha;
S = @(omega) A2*(omega.^2 + r^2).^(-p);

figure
plot(x,S(x)), xlog, ylog, hold on
scatter(j,S(j))

plot(x,H(x/j_star)/j_star)
scatter(j,H(j/j_star)/j_star)

trapz(x,S(x))
sum(S(j))
sum(H(j/j_star)/j_star)