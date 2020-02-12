clear all
clc

%% Initialize some symbolic variables
syms U V Phi M g0 g1 g2 g3 g4 g5 g6 g7 g8 omega_phase...
     k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star...
     k1g_star k2g_star k3g_star k4g_star k5g_star k6g_star k7g_star k8g_star...
     Press omega f0 f1 f2 f3 f4 f5 f6 f7 f8 Fx Fy R real

%% Define lattice directions, weights and other useful quantities of the D2Q9 model
np = 9; % number of lattice directions
cx = [0 1 0 -1 0 1 -1 -1 1];
cy = [0 0 1 0 -1 1 1 -1 -1];
w = [4/9 1/9 1/9 1/9 1/9 1/36 1/36 1/36 1/36];
cs = 1/sqrt(3);
cs2 = cs^2;
cs3 = cs^3;
cs4 = cs^4;
cs6 = cs^6;
cs8 = cs^8;

f = [f0 f1 f2 f3 f4 f5 f6 f7 f8]'; % hydrodynamic populations
feq = sym(zeros(np,1));
g = [g0 g1 g2 g3 g4 g5 g6 g7 g8]'; % phase populations
geq = sym(zeros(np,1));

Force = sym(zeros(np,1)); % generic force vector
T = sym(zeros(np,np)); %transformation matrix for central moments
M = zeros(np,np);  %transformation matrix for raw moments
K_phase = diag([1, omega_phase, omega_phase, 1, 1, 1, 1, 1, 1]); % relaxation matrices
K = diag([1, 1, 1, 1, omega, omega, 1, 1, 1]);
Id = eye(np,np); % identity matrix
for i=1:np
    % build the complete equilibria
    first_order = 1/cs2*(U*cx(i)+V*cy(i));
    second_order = 1/(2*cs4)*((cx(i)*cx(i)-1/3)*U^2+...
                              (cy(i)*cy(i)-1/3)*V^2+...
                              2*cx(i)*cy(i)*U*V);
    third_order = 1/(2*cs6)*((cx(i)^2-1/3)*cy(i)*U*U*V+(cy(i)^2-1/3)*cx(i)*U*V*V);
    fourth_order = 1/(4*cs8)*((cx(i)^2-1/3)*(cy(i)^2-1/3)*U*U*V*V);
    feq(i) = w(i)*(Press+first_order+second_order+third_order+fourth_order);
    geq(i) = Phi*w(i)*(1.+first_order+second_order+third_order+fourth_order);
    
    % build the complete forcing terms
    hat_cx = cx(i)/cs;
    hat_cy = cy(i)/cs;
    hat_ = [hat_cx, hat_cy];
    for a=1:2
        H1(a) = hat_(a);
        for b=1:2
            H2(a,b) = hat_(a)*hat_(b)-Id(a,b);
            for c=1:2
                hat_I = hat_(a)*Id(b,c)+hat_(b)*Id(a,c)+hat_(c)*Id(a,b);
                H3(a,b,c) = hat_(a)*hat_(b)*hat_(c)-hat_I;
                for d=1:2
                    hat_II = hat_(a)*hat_(b)*Id(c,d)+hat_(a)*hat_(c)*Id(b,d)+...
                             hat_(a)*hat_(d)*Id(b,c)+hat_(b)*hat_(c)*Id(a,d)+...
                             hat_(b)*hat_(d)*Id(a,c)+hat_(c)*hat_(d)*Id(a,b);
                    II = Id(a,b)*Id(c,d)+Id(a,c)*Id(b,d)+Id(a,d)*Id(b,c);
                    H4(a,b,c,d) = hat_(a)*hat_(b)*hat_(c)*hat_(d)-hat_II+II;
                end
            end
        end
    end
    first_order = 1/cs*(Fx*H1(1)+Fy*H1(2));                       
    second_order = 1/(2*cs2)*( (Fx*U+U*Fx)*H2(1,1) +...
                               (Fy*V+V*Fy)*H2(2,2) +...
                               (Fx*V+Fy*U)*H2(1,2) +...
                               (Fy*U+Fx*V)*H2(2,1) );
    third_order =  1/(6*cs3)*( (Fx*V*V+U*Fy*V+U*V*Fy)*H3(1,2,2) +...
                               (Fy*U*V+V*Fx*V+V*U*Fy)*H3(2,1,2) +... 
                               (Fy*V*U+V*Fy*U+V*V*Fx)*H3(2,2,1) +...
                               (Fx*U*V+U*Fx*V+U*U*Fy)*H3(1,1,2) +...
                               (Fy*U*U+V*Fx*U+V*U*Fx)*H3(2,1,1) +...
                               (Fx*V*U+U*Fy*U+U*V*Fx)*H3(1,2,1) );
    fourth_order = 1/(24*cs4)*( (Fx*U*V*V+U*Fx*V*V+U*U*Fy*V+U*U*V*Fy)*H4(1,1,2,2)+...
                                (Fx*V*U*V+U*Fy*U*V+U*V*Fx*V+U*V*U*Fy)*H4(1,2,1,2)+...
                                (Fx*V*V*U+U*Fy*V*U+U*V*Fy*U+U*V*V*Fx)*H4(1,2,2,1)+...
                                (Fy*U*U*V+V*Fx*U*V+V*U*Fx*V+V*U*U*Fy)*H4(2,1,1,2)+...
                                (Fy*U*V*U+V*Fx*V*U+V*U*Fy*U+V*U*V*Fx)*H4(2,1,2,1)+...
                                (Fy*V*U*U+V*Fy*U*U+V*V*Fx*U+V*V*U*Fx)*H4(2,2,1,1) );
    Force(i) = w(i)*(first_order + second_order + third_order + fourth_order);
    
    % build the transformation matrix T 
    CX = cx(i)-U;
    CY = cy(i)-V;
    T(1,i) = 1;
    T(2,i) = CX;
    T(3,i) = CY;
    T(4,i) = CX*CX+CY*CY;
    T(5,i) = CX*CX-CY*CY;
    T(6,i) = CX*CY;
    T(7,i) = CX*CX*CY;
    T(8,i) = CX*CY*CY;
    T(9,i) = CX*CX*CY*CY;
    
    % build the tranformation matrix M
    CX = cx(i);
    CY = cy(i);
    M(1,i) = 1;
    M(2,i) = CX;
    M(3,i) = CY;
    M(4,i) = CX*CX+CY*CY;
    M(5,i) = CX*CX-CY*CY;
    M(6,i) = CX*CY;
    M(7,i) = CX*CX*CY;
    M(8,i) = CX*CY*CY;
    M(9,i) = CX*CX*CY*CY;
end
T = simplify(T);
N = simplify(T*M^(-1)); %shift matrix

%% HYDRO
syms k0_pre k1_pre k2_pre k3_pre k4_pre k5_pre k6_pre k7_pre k8_pre real
syms k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star real
k_pre = simplify(T*f); %pre-collision central moments
k_eq = simplify(T*feq); % equilibrium central moments
k_force = simplify(T*Force); % forcing term central moments
k_pre(5) = k4_pre;
k_pre(6) = k5_pre;
k_star = simplify((Id-K)*k_pre + K*k_eq + (Id-K/2)*k_force); %post-collision central moments
ccode(k_star)
%post-collision populations
k_sym = [Press k1_star k2_star k3_star k4_star k5_star k6_star k7_star k8_star];
for i=1:np
    if(k_star(i)~=sym(0))
        k_star(i) = k_sym(i);
    end
end
f_post_collision_onestep = collect(simplify(T \ k_star), k_star);

% two-steps approach
raw_moments = simplify(N^(-1)*k_star);
ccode(raw_moments)
syms r0 r1 r2 r3 r4 r5 r6 r7 r8 real
r = [r0 r1 r2 r3 r4 r5 r6 r7 r8]'; %symbolic raw moments
f_post_collision_twosteps = collect(simplify(M\r),k_star);
f_in = collect(simplify(M\r),k_star);
ccode(f_in)

%% PHASE
k_pre_g = simplify(T*g); %pre-collision central moments
k_eq_g = simplify(T*geq); % equilibrium central moments
k_force_g = simplify(T*Force); % forcing term central moments
k_pre_g(2) = k1_pre;
k_pre_g(3) = k2_pre;
k_star_g = simplify((Id-K_phase)*k_pre_g + K_phase*k_eq_g + (Id-K_phase/2)*k_force_g ); %post-collision central moments
ccode(k_star_g)
%post-collision populations
syms k1_g_star k2_g_star k3_g_star k4_g_star k5_g_star k6_g_star k7_g_star k8_g_star real
k_sym_g = [Phi k1_g_star k2_g_star k3_g_star k4_g_star k5_g_star k6_g_star k7_g_star k8_g_star];
for i=1:np
    if(k_star_g(i)~=sym(0))
        k_star_g(i) = k_sym_g(i);
    end
end
g_post_collision_onestep = collect(simplify(T \ k_star_g), k_star_g);

% two-steps approach
raw_moments_g = simplify(N^(-1)*k_star_g);
ccode(raw_moments_g)
syms r0_g r1_g r2_g r3_g r4_g r5_g r6_g r7_g r8_g real
r_g = [r0_g r1_g r2_g r3_g r4_g r5_g r6_g r7_g r8_g]'; %symbolic raw moments
g_post_collision_twosteps = collect(simplify(M\r_g),k_star_g);
f_in = collect(simplify(M\r_g),k_star_g);
ccode(f_in)

        