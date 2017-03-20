% A MATLAB script to simulate a three-link, planar biped walker.
%
% This file is associated with the book Feedback Control of Dynamic
% Bipedal Robot Locomotion by Eric R. Westervelt, Jessy W. Grizzle,
% Christine Chevallereau, Jun-Ho Choi, and Benjamin Morris published
% by Taylor & Francis/CRC Press in 2007.
%
% Copyright (c) 2007 by Eric R. Westervelt, Jessy W. Grizzle, Christine
% Chevallereau, Jun-Ho Choi, and Benjamin Morris.  This code may be
% freely used for noncommercial ends.  If use of this code in part or in
% whole results in publication, proper citation must be included in that
% publication.  This code comes with no guarantees or support.
%
% Eric Westervelt
% 20 February 2007

function varargout = walker_main(t,x,flag,opt_param)

if nargin == 0
    flag = 'demo';
end

switch flag
    case ''                                 % Return dx/dt = dynamics(t,x).
        varargout{1} = f(t,x,opt_param);
    case 'events'                           % Return [value,isterminal,direction].
        [varargout{1:3}] = events(t,x);
    case 'demo'                             % Run a demo.
        varargout{1} = demo(x, opt_param, t);
    otherwise
        error(['Unknown flag ''' flag '''.']);
end


%% --------------------------------------------------------------------------
%% This is the system dynamics function

function dx = f(t,x,a)

global t_2 torque y force

[D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,dLfH] = dynamics_three_link(x(1:6),a);
Fx = inv(D)*(-C*x(4:6)-G);
Gx = inv(D)*B;

% Bernstein-Bhat controller (uses feedback linearization)
v = control_three_link(H,LfH);

% Used for controller that use feedback linearization
u = inv(dLfH*[zeros(3,2);Gx])*(v-dLfH*[x(4:6);Fx]);

dx(1:3) = x(4:6);
dx(4:6) = Fx+Gx*u;
dx = dx';

torque = [torque ; u.'];
t_2 = [t_2 ; t];
y = [y ; H.'];
[f_tan,f_norm] = stance_force_three_link(x(1:6),dx(1:6),u);
force = [force ; f_tan f_norm];


%% --------------------------------------------------------------------------
%% Locate the time when critical angle of stance leg minus stance leg angle
%% passes through zero in a decreasing direction and stop integration.

function [value,isterminal,direction] = events(t,x)

persistent control_call_cnt
if isempty(control_call_cnt) || (t == 0)
    control_call_cnt = 0;
else
    control_call_cnt = control_call_cnt + 1;
end

%% th1d  -- is the switching angle
%
if 1
    [th3d,th1d,alpha,epsilon] = control_params_three_link;
    [r,m,Mh,Mt,L,g] = model_params_three_link;
    value(1) = th1d-x(1);         % when stance leg attains angle of th1d
    value(2) = r*cos(x(1))-0.5*r; % hips get too close to ground--kill simulation
    isterminal = [1,1];           % stop when this event occurs
    direction = [-1,-1];          % decreasing direction detection
else % no events, just watch him fall!!
    value=1;
    isterminal=1;
    direction=1;
end


%% --------------------------------------------------------------------------
%% a demo function

function ret = demo(omega_1, a, isplot)

global t_2 torque y force

torque = [];
t_2 = [];
y = [];
force = [];

tstart = 0;
tfinal = 5;

%% The optimization parameters
%
%a = [0.512 0.073 0.035 -0.819 -2.27 3.26 3.11 1.89];

%omega_1 = 1.55;
x0 = sigma_three_link(omega_1,a(1,:));
x0 = transition_three_link(x0).';
x0 = x0(1:6);

options = odeset('Events','on','Refine',4,'RelTol',10^-5,'AbsTol',10^-6);

tout = tstart;
xout = x0.';
teout = []; xeout = []; ieout = [];

disp('(impact ratio is the ratio of tangential to normal');
disp('forces of the tip of the swing leg at impact)');

for i = 1:20 % run five steps
    % Solve until the first terminal event.
    [t,x,te,xe,ie] = ode45('walker_main',[tstart tfinal],x0,options,a(mod(i,3)+1,:));
    
    % Accumulate output.  tout and xout are passed out as output arguments
    nt = length(t);
    tout = [tout; t(2:nt)];
    xout = [xout;x(2:nt,:)];
    teout = [teout; te]; % Events at tstart are never reported.
    xeout = [xeout; xe];
    ieout = [ieout; ie];
    
    % Set the new initial conditions (after impact).
    x0=transition_three_link(x(nt,:));
    
    % display some useful information to the user
    disp(['step: ',num2str(i),', impact ratio:  ',num2str(x0(7)/x0(8))])
    
    % Only positions and velocities needed as initial conditions
    x0=x0(1:6);
    
    tstart = t(nt);
    if tstart>=tfinal
        break
    end
end

disp('by Eric R. Westervelt, Jessy W. Grizzle,');
disp('Christine Chevallereau, Jun-Ho Choi, and Benjamin Morris');

% if size(xout,1) > 100
%     ret = mean(abs(xout(100:end,4)))
%     if ret > 10^2 % unstable walking
%         ret = 0;
%     end
% else
%     ret = 0;
% end
%% Draw some useful graphs
if 0
    fig_hl=figure(2);
    set(fig_hl,'Position', [200 100 400 450]);
    set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
    subplot(2,1,1)
    plot(tout,xout(:,1),'-',tout,xout(:,2),'--',tout,xout(:,3),'-.')
    legend('\theta_1','\theta_2','\theta_3',-1)
    grid
    title('Joint Positions')
    subplot(2,1,2)
    plot(tout,xout(:,4),'-',tout,xout(:,5),'--',tout,xout(:,6),'-.')
    legend('\theta_1 (dot)','\theta_2 (dot)','\theta_3 (dot)',-1);
    xlabel('time (sec)')
    grid
    title('Joint Velocities')
    
    fig_hl=figure(3);
    set(fig_hl,'Position', [220 120 400 450])
    set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
    subplot(2,1,1)
    plot(t_2,force(:,1),'-b',t_2,force(:,2),'--r')
    legend('F_{tan} (N)','F_{norm} (N)',-1)
    grid
    title('Forces on End of Stance Leg')
    subplot(2,1,2)
    plot(t_2,force(:,1)./force(:,2))
    ylabel('F_{tan}/F_{norm}');
    xlabel('time (sec)')
    grid
    
    fig_hl=figure(4);
    set(fig_hl,'Position', [240 140 400 450]);
    set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
    subplot(2,1,1)
    plot(t_2,torque(:,1))
    ylabel('\tau_1 (Nm)')
    grid
    title('Control Signals')
    subplot(2,1,2)
    plot(t_2,torque(:,2))
    ylabel('\tau_2 (Nm)');
    xlabel('time (sec)')
    grid
    
    fig_hl=figure(5);
    set(fig_hl,'Position', [260 160 400 450]);
    set(fig_hl,'PaperPosition',[1.25 1.5 6 8])
    subplot(2,1,1)
    plot(t_2,y(:,1))
    ylabel('y_1')
    title('Outputs')
    grid
    subplot(2,1,2)
    plot(t_2,y(:,2))
    ylabel('y_2');
    xlabel('time (sec)')
    grid
end
% Run the animation
ret = anim(tout,xout,1/30,1,isplot);



%% --------------------------------------------------------------------------
%% the animation function

function ret = anim(t,x,ts,speed,isplot)

% gui off
%isplot = 1;

[n,m]=size(x);
[vV,vH]=hip_vel(x); % convert angles to horizontal position of hips
pH_horiz = zeros(n,1);

% Estimate hip horizontal position by estimating integral of hip velocity
for j=2:n
    pH_horiz(j)=pH_horiz(j-1)+(t(j)-t(j-1))*vH(j-1,1);
end

[te,pH_horiz]=even_sample(t,pH_horiz,1/ts);
[te,xe]=even_sample(t,x,1/ts);
[n,m]=size(xe);

q=xe(n,1:3);
[pFoot1,pFoot2,pH,pT]=limb_position(q,pH_horiz(n));
ret = pH(1) / te(n);
if te(n) < 5
    ret = pH(1)/5;
end

k=0;

q=xe(1,1:3);
[pFoot1,pFoot2,pH,pT]=limb_position(q,pH_horiz(1));

if isplot
    fig_hl=figure(1);
    set(fig_hl,'Position', [280 180 400 350]);
    set(fig_hl,'PaperPosition',[0 0 6 5]);
    clf
    anim_axis=[-2.2 2.2 -2.2 2.2];
    axis off
    axis(anim_axis)
    grid
end
anim_axis=[-2.2 2.2 -2.2 2.2];
axis off
axis(anim_axis)
% Use actual relations between masses in animation
[r,m,Mh,Mt,L,g]=model_params_three_link;
scl=0.04; % factor to scale masses
mr_legs=m^(1/3)*scl; % radius of mass for legs
mr_torso=Mt^(1/3)*scl; % radius of mass for torso
leg1_color='g';
leg2_color='r';
torso_color='b';
ground_color='k'; % a.k.a. black

% Approximate circular shape of mass
param=linspace(0,2*pi+2*pi/50,50);
xmass_legs=mr_legs*cos(param);
ymass_legs=mr_legs*sin(param);

xmass_torso=mr_torso*cos(param);
ymass_torso=mr_torso*sin(param);

% Draw ground
if isplot
    buffer=5;
    ground=line([-buffer pH_horiz(n)+buffer],[0 0]);
    set(ground,'Color',ground_color,'LineWidth',2);
    for k=-buffer:floor(pH_horiz(n)+buffer)
        ref_tick(k+buffer+1)=line([k k],[-0.1 0]);
        set(ref_tick(k+buffer+1),'Color',ground_color);
        ref_label(k+buffer+1)=text(-0.03+k,-0.2,num2str(k));
    end
    
    % Draw leg one
    leg1=line([pFoot1(1) pH(1)],[pFoot1(2) pH(2)]);
    mass1=patch(xmass_legs+(pH(1)-pFoot1(1))/2,...
        ymass_legs+(pH(2)-pFoot1(2))/2,leg1_color);
    set(mass1,'EdgeColor',leg1_color)
    set(leg1,'LineWidth',2,'Color',leg1_color);
    
    % Draw leg two
    leg2=line([pFoot2(1) pH(1)],[pFoot2(2) pH(2)]);
    mass2=patch(xmass_legs+pH(1)-(pH(1)-pFoot2(1))/2,...
        ymass_legs+pH(2)-(pH(2)-pFoot2(2))/2,leg2_color);
    set(mass2,'EdgeColor',leg2_color)
    set(leg2,'LineWidth',2,'Color',leg2_color);
    
    % Draw torso
    torso=line([pH(1) pT(1)],[pH(2)*2 pT(2)*2]);
    torso_mass=patch(xmass_torso+pT(1),ymass_torso+pT(2),torso_color);
    set(torso_mass,'EdgeColor',torso_color)
    set(torso,'LineWidth',2,'Color',torso_color);

    for k=2:n
        q=xe(k,1:3);
        [pFoot1,pFoot2,pH,pT]=limb_position(q,pH_horiz(k));
        ret = pH(1) / te(k);
        if te(k) < 5
            ret = pH(1)/5;
        end
        
        set(leg1,'XData',[pFoot1(1) pH(1)],'YData',[pFoot1(2) pH(2)]);
        
        set(mass1,'XData',xmass_legs+(pH(1)-pFoot1(1))/2+pH_horiz(k),...
            'YData',ymass_legs+(pH(2)-pFoot1(2))/2);
        
        set(leg2,'XData',[pFoot2(1) pH(1)],'YData',[pFoot2(2) pH(2)]);
        
        set(mass2,'XData',xmass_legs+pH(1)-(pH(1)-pFoot2(1))/2,...
            'YData',ymass_legs+pH(2)-(pH(2)-pFoot2(2))/2);
        
        set(torso,'XData',[pH(1) pT(1)+(pT(1)-pH(1))],...
            'YData',[pH(2) pT(2)+(pT(2)-pH(2))]);
        
        set(torso_mass,'XData',xmass_torso+pT(1),'YData',ymass_torso+pT(2));
        
        title(['T_{est} = ',num2str(te(k),'%.1f')])
        
        
        
        new_axis=anim_axis+[pH(1) pH(1) 0 0];
        axis(new_axis);
        for j=1:length(ref_label)
            if (j-buffer-1.05<new_axis(1)) || (j-buffer-1>new_axis(2))
                set(ref_label(j),'Visible','off')
                set(ref_tick(j),'Visible','off')
            else
                set(ref_label(j),'Visible','on');
                set(ref_tick(j),'Visible','on')
            end
        end
        
        drawnow;
        pause(ts*speed);
    end
end


%% --------------------------------------------------------------------------
%% a function to calculate hip velocity

function [vV,vH] = hip_vel(x)

vV=zeros(length(x),1);
vH=cos(x(:,1)).*x(:,4); % estimate of horizontal velocity of hips


%% --------------------------------------------------------------------------
%% a function to calculate the limb position

function [pFoot1,pFoot2,pH,pT] = limb_position(q,pH_horiz)

% Use position of hips as location of stance leg foot.
[r,m,Mh,Mt,L,g]=model_params_three_link;

pFoot1=[pH_horiz; 0];
pH=[pFoot1(1)+r*sin(q(1)); pFoot1(2)+r*cos(q(1))];
pFoot2=[pH(1)-r*sin(q(2)); pH(2)-r*cos(q(2))];
pT=[pH(1)+L*sin(q(3)); pH(2)+L*cos(q(3))];


%% --------------------------------------------------------------------------
%% CONVERTS A RANDOMLY SAMPLED SIGNAL SET INTO AN EVENLY SAMPLED SIGNAL SET
%% (by interpolation)
%%
%% written by Haldun Komsuoglu, 7/23/1999

function [Et, Ex] = even_sample(t, x, Fs)

% Obtain the process related parameters
N = size(x, 2);    % number of signals to be interpolated
M = size(t, 1);    % Number of samples provided
t0 = t(1,1);       % Initial time
tf = t(M,1);       % Final time
EM = (tf-t0)*Fs;   % Number of samples in the evenly sampled case with
% the specified sampling frequency
Et = linspace(t0, tf, EM)';

% Using linear interpolation (used to be cubic spline interpolation)
% and re-sample each signal to obtain the evenly sampled forms
for s = 1:N,
    Ex(:,s) = interp1(t(:,1), x(:,s), Et(:,1));
end
