clear, clc

a = 0.1;
b = 0.2;
omega = -1;
eps = 1e-8;
t = linspace(0, 1, 101);
pos = zeros(2, length(t));
vel = zeros(2, length(t));

%initial guesses
x_init = [0.5; 0.5];
v_init = [0.5; 0.5];

%position analysis

for ii = 1:length(t)
    
    fi = pi/6 + omega * t(ii);
    
    f =  @(x) [a*cos(fi) + b*cos(x(1)) - x(2);
               a*sin(fi) - b*sin(x(1))];
           
    Jacobi = @(x) [-b*sin(x(1)), -1;
                   -b*cos(x(1)), 0];

    [pos(:,ii), iter] = NR_method(f, Jacobi, x_init, eps);

    theta = pos(1,:);
    d = pos(2,:); 
end

    
figure
set(gcf,'position',[200,200,1000,300])
subplot(1,2,1)
plot(t, theta, 'LineWidth', 1.5)
xlabel('t [s]', 'Interpreter','latex')
ylabel('\theta [rad]')
grid on
axis tight

subplot(1,2,2)
plot(t, d, 'r','LineWidth', 1.5)
xlabel('t [s]', 'Interpreter','latex')
ylabel('d [m]', 'Interpreter','latex')
grid on
axis tight

%velocity analysis

for ii = 1:length(t)
    
    fi = pi/6 + omega * t(ii);
    fi_dot = omega;

    f_dot = @(x) [-a*sin(fi)*fi_dot - b*sin(theta(ii))*x(1) - x(2);
                    a*cos(fi)*fi_dot - b*cos(theta(ii))*x(1)];

    Jacobi_dot = @(x) [-b*sin(theta(ii)), -1;
                  -b*cos(theta(ii)), 0];

    [vel(:,ii), iter2] = NR_method(f_dot, Jacobi_dot, v_init, eps);

    theta_dot = vel(1,:);
    d_dot = vel(2,:);
end

figure
set(gcf,'position',[200,200,1000,300])
subplot(1,2,1)
plot(t, theta_dot, 'LineWidth', 1.5)
xlabel('t [s]', 'Interpreter','latex')
ylabel('$\dot{\theta}$ [rad/s]', 'Interpreter', 'latex')
grid on
axis tight

subplot(1,2,2)
plot(t, d_dot, 'r','LineWidth', 1.5)
xlabel('t [s]', 'Interpreter','latex')
ylabel('$\dot{d}$ [m/s]', 'Interpreter', 'latex')
grid on
axis tight