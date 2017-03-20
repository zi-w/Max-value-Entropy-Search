%gen_model_three_link.m
%
%   Generate model for three-link biped walker and output to
%   an MATLAB funtion file dynamics.m
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

clear

filename = strrep(mfilename,'gen_model','model_save');

if isempty(dir([filename,'.mat']))
	%if 1
	clear

	filename = strrep(mfilename,'gen_model','model_save');

	%% variable declarations
	%
	syms th1 th2 th3 z1 z2 real
	syms dth1 dth2 dth3 dz1 dz2 real
	syms m Mh Mt g L r real
	syms th1d th3d
	syms a01 a11 a21 a31
	syms a02 a12 a22 a32

	%% generalized coordinates
	%
	q=[th1;th2;th3];
	qe=[q;z1;z2];

	%% first derivative of generalized coordinates
	%
	dq=[dth1;dth2;dth3];
	dqe=[dq;dz1;dz2];

	%% Generate De matrix.  To do so, we must go through the trouble of
	%% defining the positions of the masses using the augmented
	%% configuration variables.
	%

	%% position of masses in system
	%
	p_m1=[z1+r/2*sin(th1); z2+r/2*cos(th1)];
	p_Mh=[z1+r*sin(th1); z2+r*cos(th1)];
	p_Mt=p_Mh+[L*sin(th3); L*cos(th3)];
	p_m2=p_Mh-[r/2*sin(th2); r/2*cos(th2)];

	%% velocities of masses in system
	%
	v_m1=jacobian(p_m1,qe)*dqe;
	v_Mh=jacobian(p_Mh,qe)*dqe;
	v_Mt=jacobian(p_Mt,qe)*dqe;
	v_m2=jacobian(p_m2,qe)*dqe;

	%% kinetic energy of masses in system
	%
	KE_m1=simple(m/2*v_m1'*v_m1);
	KE_Mh=simple(Mh/2*v_Mh'*v_Mh);
	KE_Mt=simple(Mt/2*v_Mt'*v_Mt);
	KE_m2=simple(m/2*v_m2'*v_m2);

	%% total kinetic energy of system
	%
	KE=KE_m1+KE_Mh+KE_Mt+KE_m2;

	%% potential energy of masses in system
	%
	PE_m1=p_m1(2)*m*g;
	PE_Mh=p_Mh(2)*Mh*g;
	PE_Mt=p_Mt(2)*Mt*g;
	PE_m2=p_m2(2)*m*g;

	%% total potential energy of system
	%
	PE=PE_m1+PE_Mh+PE_Mt+PE_m2;

	De=simplify(jacobian(jacobian(KE,dqe).',dqe));

	N=max(size(qe));
    syms Ce
	for k=1:N,
		for j=1:N,
			Ce(k,j)=0*g;
			for i=1:N,
				Ce(k,j)=Ce(k,j)+1/2*(diff(De(k,j),qe(i))+...
					diff(De(k,i),qe(j))-...
					diff(De(i,j),qe(k)))*dqe(i);
			end
		end
	end

	Ge=jacobian(PE,qe).';

	%% Now that We're done generating the De matrix we now can use the
	%% non-augmented configuration variables to derive the required
	%% matrices.
	%

	%% position of masses in system
	%
	p_m1=[r/2*sin(th1); r/2*cos(th1)];
	p_Mh=[r*sin(th1); r*cos(th1)];
	p_Mt=p_Mh+[L*sin(th3); L*cos(th3)];
	p_m2=p_Mh-[r/2*sin(th2); r/2*cos(th2)];

	%% velocities of masses in system
	%
	v_m1=jacobian(p_m1,q)*dq;
	v_Mh=jacobian(p_Mh,q)*dq;
	v_Mt=jacobian(p_Mt,q)*dq;
	v_m2=jacobian(p_m2,q)*dq;

	%% kinetic energy of masses in system
	%
	KE_m1=simple(m/2*v_m1'*v_m1);
	KE_Mh=simple(Mh/2*v_Mh'*v_Mh);
	KE_Mt=simple(Mt/2*v_Mt'*v_Mt);
	KE_m2=simple(m/2*v_m2'*v_m2);

	%% total kinetic energy of system
	%
	KE=KE_m1+KE_Mh+KE_Mt+KE_m2;

	%% potential energy of masses in system
	%
	PE_m1=p_m1(2)*m*g;
	PE_Mh=p_Mh(2)*Mh*g;
	PE_Mt=p_Mt(2)*Mt*g;
	PE_m2=p_m2(2)*m*g;

	%% total potential energy of system
	%
	PE=PE_m1+PE_Mh+PE_Mt+PE_m2;

	%% the Lagragian
	%
	%Lag=KE-PE;

	%% Form D, C, G, B, and F matrices of
	%%
	%% D(q)ddq+C(q,dq)dq+G(q)=B*tau
	%%
	%% where tau=[tau1; tau2]
	%
	D=jacobian(jacobian(KE,dq).',dq);

	N=max(size(q));
    syms C
	for k=1:N,
		for j=1:N,
			C(k,j)=0*g;
			for i=1:N,
				C(k,j)=C(k,j)+1/2*(diff(D(k,j),q(i))+...
					diff(D(k,i),q(j))-...
					diff(D(i,j),q(k)))*dq(i);
			end
		end
	end

	G=jacobian(PE,q).';

	%% change of coordinates for torque inputs
	%
	T=[-1  0  1;
		0 -1  1;
		0  0  1]*q;

	B=jacobian(T,q)'*[1 0;
		0 1;
		0 0];

	Psi=[z1 + r*sin(th1)-r*sin(th2);
		z2 + r*cos(th1)-r*cos(th2)];

	%% Used to calculate post-impact conditions
	%
	E=jacobian(Psi,qe);

	%% output objective
	%
	H=[th3-th3d; th1+th2];
	LfH=jacobian(H,q)*dq;  % there is a simplification because H does
	% not involve derivates of the gen. coords.
	dLfH=jacobian(LfH,[q,dq]); % need to multiply this quantity times
	% g(x)=[zeros(3,2);inv(D)*B]; to get
	% LgLfH, and by
	% f(x)=[eye(3)*dq;inv(D)*(-C*dq-G)]

	%% output function for use in optimization
	%
	Ha=[th3-(a01+a11*th1+a21*th1^2+a31*th1^3);
		th2-(-th1+(a02+a12*th1+a22*th1^2+a32*th1^3)*(th1-th1d)*(th1+th1d))];
	LfHa=jacobian(Ha,q)*dq;
	dLfHa=jacobian(LfHa,[q,dq]); % need to multiply this quantity times
	% g(x)=[zeros(3,2);inv(D)*B]; to get
	% LgLfH, and by
	% f(x)=[eye(3)*dq;inv(D)*(-C*dq-G)]

	%% Feedback gain for linearized system
	%
	Al=[0 0 1 0;
		0 0 0 1;
		0 0 0 0;
		0 0 0 0];
	Bl=[0 0;
		0 0;
		1 0;
		0 1];
	K=place(Al,Bl,[-10;-11;-12;-13]);

	P=lyap(Al-Bl*K,diag([10,10,10,10]));

	V=[H;LfH]'*P*[H;LfH];

	dV=jacobian(V,[q dq]);
	dVl=2*[H;LfH]'*P;

	save(filename);
else
	clear
	load(strrep(mfilename,'gen_model','model_save'));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Output functions for use in simulation
%%
%

N=max(size(q));

%% First, output model to a file called
%%
%%    dynamics+<model type>.m
%%
%

%% Output file header
%
fcn_name=strrep(mfilename,'gen_model','dynamics');
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,dLfH]=%s',fcn_name);
fprintf(fid,'(x,a)\n');
fprintf(fid,'%% %s    Model of three-link biped walker model.\n',...
	upper(fcn_name));
fprintf(fid,'%%    [D,C,G,B,K,dV,dVl,Al,Bl,H,LfH,DLFH] = %s',upper(fcn_name));
fprintf(fid,'(X,\n%%      A) is the three-link\n');
fprintf(fid,'%%    biped walking model. (x is of dimension %s)\n\n',...
	num2str(2*N));

fprintf(fid,'%% Eric Westervelt\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[r,m,Mh,Mt,L,g]=model_params_three_link;\n\n');
fprintf(fid,'[th3d,th1d,alpha,epsilon]=control_params_three_link;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'th1=x(1); th2=x(2); th3=x(3);\n');
fprintf(fid,'dth1=x(4); dth2=x(5); dth3=x(6);\n\n');

%% Model output
%
fprintf(fid,'%% D matrix\n');
fprintf(fid,'D=zeros(%s);\n',num2str(N));
for k=1:N,
	for j=1:N,
		if D(k,j)~=0
			ttt=char(D(k,j));
			fprintf(fid,'D(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% C matrix\n');
fprintf(fid,'C=zeros(%s);\n',num2str(N));
for k=1:N,
	for j=1:N,
		if C(k,j)~=0
			ttt=char(C(k,j));
			fprintf(fid,'C(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% G matrix\n');
fprintf(fid,'G=zeros(%s,1);\n',num2str(N));
for k=1:N,
	if G(k)~=0
		ttt=char(G(k));
		fprintf(fid,'G(%s)=%s;\n',num2str(k),ttt);
	end
end

fprintf(fid,'\n%% B matrix\n');
[N,M]=size(B);
fprintf(fid,'B=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if B(k,j)~=0
			ttt=char(B(k,j));
			fprintf(fid,'B(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

%% Temporary stuff
%
fprintf(fid,'\n%% K matrix\n');
[N,M]=size(K);
fprintf(fid,'K=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if K(k,j)~=0
			fprintf(fid,'K(%s,%s)=%s;\n',num2str(k),num2str(j),num2str(K(k,j)));
		end
	end
end

fprintf(fid,'\n%% dV matrix\n');
[N,M]=size(dV);
fprintf(fid,'dV=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if dV(k,j)~=0
			ttt=char(dV(k,j));
			fprintf(fid,'dV(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% dVl matrix\n');
[N,M]=size(dVl);
fprintf(fid,'dVl=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if dVl(k,j)~=0
			ttt=char(dVl(k,j));
			fprintf(fid,'dVl(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% Al matrix\n');
[N,M]=size(Al);
fprintf(fid,'Al=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if Al(k,j)~=0
			fprintf(fid,'Al(%s,%s)=%s;\n',num2str(k),num2str(j),num2str(Al(k,j)));
		end
	end
end

fprintf(fid,'\n%% Bl matrix\n');
[N,M]=size(Bl);
fprintf(fid,'Bl=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if Bl(k,j)~=0
			fprintf(fid,'Bl(%s,%s)=%s;\n',num2str(k),num2str(j),num2str(Bl(k,j)));
		end
	end
end

%% Reassign optimization parameters
%
fprintf(fid,'a01=a(1); a11=a(2); a21=a(3); a31=a(4);\n');
fprintf(fid,'a02=a(5); a12=a(6); a22=a(7); a32=a(8);\n\n');

fprintf(fid,'%% Ha matrix\n');
[N,M]=size(Ha);
fprintf(fid,'H=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if Ha(k,j)~=0
			ttt=char(Ha(k,j));
			fprintf(fid,'H(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end
fprintf(fid,'\n%% LfHa matrix\n');
[N,M]=size(LfHa);
fprintf(fid,'LfH=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if LfHa(k,j)~=0
			ttt=char(LfHa(k,j));
			fprintf(fid,'LfH(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end
fprintf(fid,'\n%% dLfHa matrix\n');
[N,M]=size(dLfHa);
fprintf(fid,'dLfH=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if dLfHa(k,j)~=0
			ttt=char(dLfHa(k,j));
			fprintf(fid,'dLfH(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end
fclose(fid);


%% Second, output model to a file called
%%
%%    transition+<model type>.m
%%
%

N=max(size(q));

%% Output file header
%
fcn_name=strrep(mfilename,'gen_model','transition');
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [x_new,z2_new]=%s(x)\n',fcn_name);
fprintf(fid,'%% %s    Calculate the state of the system after impact.\n',...
	upper(fcn_name));
fprintf(fid,'%%          (Last two entries are the forces at the toe.\n');
fprintf(fid,'%%    [X_NEW,Z2_NEW] = %s(X) is the transition function for\n',...
	upper(fcn_name));
fprintf(fid,'%%    biped walking model. (x is of dimension %s)\n\n',...
	num2str(2*N+2));
fprintf(fid,'%% Eric Westervelt\n');
fprintf(fid,'%% %s\n\n',datestr(now));

N=max(size(qe));

%% Read in constants
%
fprintf(fid,'[r,m,Mh,Mt,L,g]=model_params_three_link;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'th1=x(1); th2=x(2); th3=x(3);\n\n');

fprintf(fid,'%% De matrix\n');
fprintf(fid,'De=zeros(%s,%s);\n',num2str(N),num2str(N));
for k=1:N,
	for j=1:N,
		if De(k,j)~=0
			ttt=char(De(k,j));
			fprintf(fid,'De(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

[N,M]=size(E);
fprintf(fid,'\n%% E matrix\n');
fprintf(fid,'E=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if E(k,j)~=0
			ttt=char(E(k,j));
			fprintf(fid,'E(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% See Grizzle''s paper, page 28 for equation...\n');
fprintf(fid,...
	'tmp_vec=inv([De -E'';E zeros(2)])*[De*[x(4:6)'';zeros(2,1)];zeros(2,1)];\n\n');

fprintf(fid,'x_new(1)=x(2);\n');
fprintf(fid,'x_new(2)=x(1);\n');
fprintf(fid,'x_new(3)=x(3);\n');
fprintf(fid,'x_new(4)=tmp_vec(2);\n');
fprintf(fid,'x_new(5)=tmp_vec(1);\n');
fprintf(fid,'x_new(6)=tmp_vec(3);\n');
fprintf(fid,'x_new(7)=tmp_vec(6);\n');
fprintf(fid,'x_new(8)=tmp_vec(7);\n');
fprintf(fid,'z2_new=tmp_vec(5);\n');
fclose(fid);


%% Third, output model to a file called
%%
%%    control+<model type>.m
%%
%

%% Output file header
%
fcn_name=strrep(mfilename,'gen_model','control');
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [v]=%s(H,LfH)\n',fcn_name);
fprintf(fid,'%% %s    Calculate the control.\n',...
	upper(fcn_name));
fprintf(fid,'%%    [V] = %s(X,H,LFH) is the control for the\n',...
	upper(fcn_name));
fprintf(fid,'%%    feedback linearized biped walking model.\n');
fprintf(fid,'%%\n\n');
fprintf(fid,'%% Eric Westervelt\n');
fprintf(fid,'%% %s\n\n',datestr(now));

%% Read in constants
%
fprintf(fid,'[th3d,th1d,alpha,epsilon]=control_params_three_link;\n\n');

fprintf(fid,'%% LfH scaling\n');
fprintf(fid,'LfH=epsilon*LfH;\n\n');

fprintf(fid,'%% phi fcns\n');
fprintf(fid,'phi1=H(1)+1/(2-alpha)*sign(LfH(1))*abs(LfH(1))^(2-alpha);\n');
fprintf(fid,'phi2=H(2)+1/(2-alpha)*sign(LfH(2))*abs(LfH(2))^(2-alpha);\n\n');

fprintf(fid,'%% psi fcns\n');
fprintf(fid,'psi(1,1)=-sign(LfH(1))*abs(LfH(1))^alpha...\n');
fprintf(fid,'         -sign(phi1)*abs(phi1)^(alpha/(2-alpha));\n');
fprintf(fid,'psi(2,1)=-sign(LfH(2))*abs(LfH(2))^alpha...\n');
fprintf(fid,'         -sign(phi2)*abs(phi2)^(alpha/(2-alpha));\n\n');

fprintf(fid,'%% calculate control\n');
fprintf(fid,'v=1/epsilon^2*psi;\n');
fclose(fid);


%% Fourth, stance-leg force calculation routine (for impacts) to a file called
%%
%%    stance_force_+<model type>.m
%%
%

N=max(size(q));
Ne=max(size(qe));

%% Output file header
%
fcn_name=strrep(mfilename,'gen_model','stance_force');
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [f_tan,f_norm]=%s(x,dx,u)\n',fcn_name);
fprintf(fid,'%% %s    Calculate the forces on the stance\n',...
	upper(fcn_name));
fprintf(fid,'%%                            leg during impact.\n');
fprintf(fid,'%%    [F_TAN,F_NORM] = %s(X,DX,U) are the forces on the\n',...
	upper(fcn_name));
fprintf(fid,'%%    stance leg at impact.\n\n');
fprintf(fid,'%% Eric Westervelt\n');
fprintf(fid,'%% %s\n\n',datestr(now));

% Read in constants
%
fprintf(fid,'[r,m,Mh,Mt,L,g]=model_params_three_link;\n\n');
fprintf(fid,'[th3d,th1d,alpha,epsilon]=control_params_three_link;\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'th1=x(1); th2=x(2); th3=x(3);\n');
fprintf(fid,'dth1=x(4); dth2=x(5); dth3=x(6);\n\n');

%% Model output
%
fprintf(fid,'%% De11 matrix\n');
fprintf(fid,'De11=zeros(%s,%s);\n',num2str(N),num2str(N));
for k=1:N,
	for j=1:N,
		if De(k,j)~=0
			ttt=char(De(k,j));
			fprintf(fid,'De11(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% De12 matrix\n');
fprintf(fid,'De12=zeros(%s,%s);\n',num2str(N),num2str(Ne-N));
for k=1:N,
	for j=1:Ne-N,
		if De(k,j+N)~=0
			ttt=char(De(k,j+N));
			fprintf(fid,'De12(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% De22 matrix\n');
fprintf(fid,'De22=zeros(%s,%s);\n',num2str(Ne-N),num2str(Ne-N));
for k=1:Ne-N,
	for j=1:Ne-N,
		if De(k+N,j+N)~=0
			ttt=char(De(k+N,j+N));
			fprintf(fid,'De22(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% Ce11 matrix\n');
fprintf(fid,'Ce11=zeros(%s,%s);\n',num2str(N),num2str(N));
for k=1:N,
	for j=1:N,
		if Ce(k,j)~=0
			ttt=char(Ce(k,j));
			fprintf(fid,'Ce11(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% Ce21 matrix\n');
fprintf(fid,'Ce21=zeros(%s,%s);\n',num2str(Ne-N),num2str(N));
for k=1:Ne-N,
	for j=1:N,
		if Ce(k+N,j)~=0
			ttt=char(Ce(k+N,j));
			fprintf(fid,'Ce21(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% Ge1 matrix\n');
fprintf(fid,'Ge1=zeros(%s,%s);\n',num2str(N),num2str(1));
for k=1:N,
	if Ge(k)~=0
		ttt=char(Ge(k));
		fprintf(fid,'Ge1(%s,%s)=%s;\n',num2str(k),num2str(1),ttt);
	end
end

fprintf(fid,'\n%% Ge2 matrix\n');
fprintf(fid,'Ge2=zeros(%s,%s);\n',num2str(Ne-N),num2str(1));
for k=1:Ne-N,
	if Ge(k+N)~=0
		ttt=char(Ge(k+N));
		fprintf(fid,'Ge2(%s,%s)=%s;\n',num2str(k),num2str(1),ttt);
	end
end

fprintf(fid,'\n%% B matrix\n');
[N,M]=size(B);
fprintf(fid,'B=zeros(%s,%s);\n',num2str(N),num2str(M));
for k=1:N,
	for j=1:M,
		if B(k,j)~=0
			ttt=char(B(k,j));
			fprintf(fid,'B(%s,%s)=%s;\n',num2str(k),num2str(j),ttt);
		end
	end
end

fprintf(fid,'\n%% See my notes, 2/16/200 for equations...\n');

fprintf(fid,'DD=inv((De12*inv(De22)).''*De12*inv(De22))*(De12*inv(De22)).'';\n');

fprintf(fid,'F=DD*(-(De11-De12*inv(De22)*De12.'')...\n');
fprintf(fid,'  *dx(4:6)+(De12*inv(De22)*Ce21-Ce11)...\n');
fprintf(fid,'  *dx(1:3)+De12*inv(De22)*Ge2-Ge1+B*u);\n\n');

fprintf(fid,'f_tan=F(1);\n');
fprintf(fid,'f_norm=F(2);\n');
fclose(fid);


%% Fifth, sigma function that maps theta to state of biped just before impact.
%%
%%    sigma_+<model type>.m
%%
%

N=max(size(q));
Ne=max(size(qe));

%% Output file header
%
fcn_name=strrep(mfilename,'gen_model','sigma');
fid=fopen([fcn_name,'.m'],'w');
fprintf(fid,'function [x]=%s(omega_1_minus,a)\n',fcn_name);
fprintf(fid,'%% %s    Maps velocity of stance leg just before\n',...
	upper(fcn_name));
fprintf(fid,'%%         impact to state of the system just before impact.\n');
fprintf(fid,'%%    [X] = %s(OMEGA_1_MINUS,A)\n\n',upper(fcn_name));
fprintf(fid,'%% Eric Westervelt\n');
fprintf(fid,'%% %s\n\n',datestr(now));

% Read in constants
%
fprintf(fid,'[th3d,th1d,alpha,epsilon]=control_params_three_link;\n\n');

%% Reassign optimization parameters
%
fprintf(fid,'a01=a(1); a11=a(2); a21=a(3); a31=a(4);\n');
fprintf(fid,'a02=a(5); a12=a(6); a22=a(7); a32=a(8);\n\n');

%% Reassign configuration parameters
%
fprintf(fid,'th1=th1d;\n');
fprintf(fid,'dth1=omega_1_minus;\n\n');

syms th3 dth2 dth3

ttt=char(solve(Ha(1),th3));
fprintf(fid,'th3 = %s;\n',ttt);
ttt=char(solve(LfHa(2),dth2));
fprintf(fid,'dth2 = %s;\n',ttt);
ttt=char(solve(LfHa(1),dth3));
fprintf(fid,'dth3 = %s;\n\n',ttt);

fprintf(fid,'x = [th1,-th1,th3,dth1,dth2,dth3];');
fclose(fid);