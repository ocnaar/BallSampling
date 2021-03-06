%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octavio Narvaez-Aroche                                                  %
% ocnaar@berkeley.edu                                                     %
% Berkeley Center for Control and Identification                          %
% Fall 2017                                                               %
%                                                                         %
% Random sample of points from 2-dimensional L-p balls of different radii.%
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of samples. 
m = 500;

% Radii.
r1 = 21.5;
r2 = 17;
rinf = 10.5;

% Center 
x = [0;0];

% Sample points in L-p balls. 
L1 = nBallSampling(m,r1,x,1,1);
L2 = nBallSampling(m,r2,x,2,2);
L3 = nBallSampling(m,rinf,x,Inf,3);

% Plot samples, and boundaries of L-p balls.
figure()
plot(L1(1,:),L1(2,:),'rx',L2(1,:),L2(2,:),'bx',L3(1,:),L3(2,:),'gx','LineWidth',2)
hold on
plot([r1,0,-r1,0,r1],[0,r1,0,-r1,0],'r-',r2*cos(0:0.01:2*pi),r2*sin(0:0.01:2*pi),'b-',[rinf,-rinf,-rinf,rinf,rinf],[rinf,rinf,-rinf,-rinf,rinf],'g-','LineWidth',2)
axis equal
grid
xlabel('x_1')
ylabel('x_2')
legend({'$\left\Vert x \right\Vert_1 \leq r_1$','$\left\Vert x \right\Vert_2 \leq r_2$','$\left\Vert x \right\Vert_\infty \leq r_\infty$'},'Interpreter','Latex')