clear;

P.M = 1;
P.K = 3000;
P.R = 1;

x0 = [0; 0];
f = 1000;

[t, y] = ode45(@(t, x)rhs_fun(t, x, f, P), [0 3], x0);

%
figure;
plot(t, y(:,2));