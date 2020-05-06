function stiff_ode_example(dt, tfinal)

% This example numerically solves the stiff ode y'(t) = -15*y(t) using Euler's
% method. This has an exact solution given by y(t) = exp(-15*t). The inputs
% are dt, the time step size, and tfinal, the final simulation time.

N = ceil(tfinal/dt);    %Determine the number of time steps
y = zeros(N);           %Initialize y to hold the numerical solution
y_exact = zeros(N);     %Initialize y_exact to hold the exact solution
time = zeros(N);        %Initialize time to hold the times

y(1) = 1;               %Initial conditions
y_exact(1) = 1;
time(1) = 0;

%Loop through all time steps and calculate y using Euler's method
for i = 2:N,
  y(i)= y(i-1)+dt*(-15*y(i-1));
  time(i)=(i-1)*dt;                 %update time
  y_exact(i) = exp(-15*time(i));    %save exact solution
end

%graph the numerical and the exact solutions
hold off
plot(time, y, '-b');
hold on
plot(time, y_exact, '-r');
xlabel('time');
ylabel('y');
legend('Eulers method', 'Exact solution');