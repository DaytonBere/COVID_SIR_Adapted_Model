% Code below is to run epidemic models -- each block starting with %% is a%
% separate model.   This script requires the short function   sir_ode.m
%
%%  Basic SIR model
% $$\eqalign{
% \frac{dS}{dt} &= -\alpha\,S\,Ia -\beta\,S\,Is \cr
% \frac{dIa}{dt} &= \alpha\,S\,Ia - kIa +Is(lambda1 - lambda2)\cr
% \frac{dIs}{dt} &= \beta\,S\,Is - kIs + Ia(lambda2 - lambda1) \cr
% \frac{dR}{dt} &= k(Ia + Is) \cr
% }$$

  T = 70 ;  % final time period
  beta = 0.35;
  alpha = 0.3;
  lambda = 0.05;
  k1 = 0.15; 
  k2 = 0.14;
  y0 = [0.96 0.2 0.2 0];   % initial state for compartments
  
  f = @(t,y) [-alpha*y(1)*y(2)-beta*y(1)*y(3);
                alpha*y(1)*y(2)-k1*y(2)-y(2)*lambda;
                beta*y(1)*y(3)-k2*y(3)+y(2)*lambda;
                k1*y(2)+k2*y(3)];
                
  [t,yvals] = ode45(f,[0 T], y0);   % actual ode solver
   figure(1); plot(t,yvals(:,1),'g',t,yvals(:,2),'r',t,yvals(:,3),'b',t,yvals(:,4),'c');
   legend('S','Ia','Is','R', 'Location','northeast');  title('SIIR Model');
  
% y0 = [0.96 0.3 0.01 0];
% tspan = [0 T]; 
% [tvals,yvals] = ode45(@system1,tspan,y0);
% plot(tvals,yvals,'-o')
% 
% 
% function dydt= system1(t,y)
% dydt = zeros(4,1);
% dydt(1) = -alpha*y(1)*y(2)-beta*y(1)*y(3);
% dydt(2) = alpha*y(1)*y(2)-k*y(2)+y(3)*(lambda1-lambda2);
% dydt(3) = beta*y(1)*y(3)-k*y(3)+y(2)*(lambda2-lambda1);
% dydt(4) = k*(y(2)+y(3));
% end  



%%  SIRS model
%  $$\eqalign{
% \frac{dS}{dt} &= -\beta\,S\,I + \eta R \cr
% \frac{dI}{dt} &= \beta\,S\,I - kI \cr
% \frac{dR}{dt} &= kI - \eta R \cr
% } $$

%  T = 100 ;  % final time period
%  beta = 0.3;
%  k = 0.2; 
%  eta = 0.05;
%  y0 = [0.98 0.02 0];    % initial state for compartments 
%  f = @(t,y) [-beta*y(1)*y(2)+eta*y(3);beta*y(1)*y(2)-k*y(2);k*y(2)-eta*y(3)];
%  [t,yvals] = ode45(f,[0 T], y0);  % actual ode solver
%  figure(2);  plot(t,yvals(:,1),'g',t,yvals(:,2),'r',t,yvals(:,3),'b');
%  legend('S','I','R', 'Location','northeast');  title('SIRS model');
% %%  SEIRD model
% % $$\eqalign{
% % S' &= -\beta\,S\,I \cr
% % E' &= \beta\,S\,I - \sigma E \cr
% % I' &= \sigma\,E - (1-\alpha)k\,I - \alpha\rho\,I\cr
% % R' &= (1-\alpha)k\,I \cr
% % D' &= \alpha\rho\,I \cr
% % }$$
% 
  T = 300 ;  % final time period
  alpha = 0.6;
  beta = 0.2;
  k1 = 0.15; 
  k2 = 0.1;
  eta = 0.7;
  tau = 0.1;
  lambda = 0.05;
  sigma1 = 0.3;
  sigma2 = 0.35;
  rho = 0.07;
  y0 = [0.9 0.02 0.04 0.04 0 0];  % initial state for compartments
  f = @(t,y) [-alpha*y(1)*y(3) - beta*y(1)*y(4)+eta*y(5);
               alpha*y(1)*y(3) + beta*y(1)*y(4) - y(2)*(sigma1+sigma2);
               sigma1*y(2)-(1-tau)*k1*y(3)-tau*rho*y(3) - y(3)*lambda;
               sigma2*y(2)-(1-tau)*k2*y(4)-tau*rho*y(4) + y(3)*lambda;
               (1-tau)*(k1*y(3)+k2*y(4))-eta*y(5); 
               tau*rho*(y(3)+y(4)) ];
  [t,yvals] = ode45(f,[0 T], y0);    % actual ode solver
  figure(3);  plot(t,yvals(:,1),'c',t,yvals(:,2),'y',t,yvals(:,3),'r--');
  hold on; plot(t,yvals(:,4),'b--',t,yvals(:,5),'g',t,yvals(:,6), 'k--');
 legend('S','E','Ia', 'Is','R','D', 'Location','northeast');  title('SEIIRDS model');