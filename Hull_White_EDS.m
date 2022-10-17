tic;
% clc;
warning off;
format long;
S_0=100;%initial asset price
r=0.02;q=0;kappa=0.1;eta=2;rho=-0.7;v_0=0.05;%parameter of Hull-White model 
K=[90;100;110];% strike price
% K=100;
Tmat=1;N_monitor=4;delta_t=Tmat/N_monitor;%maturity and the number of monitor dates
V_c_benchmark=[11.0034;1.8980;0.0272];
%% grid of state variable
L=150;
n_grid=3500;
n_l=n_grid;n_r=n_grid;
delta_omega=L/n_r;
xi_m=1.1;xi_p=1.4;% dampening factor belongs to (1,1/(1-rho^2))
omega_m=(-n_l:n_r)'*delta_omega+1i*xi_m;omega_p=(-n_l:n_r)'*delta_omega+1i*xi_p;
omegam1=omega_m(n_l+1:end);
omegap1=omega_p(n_l+1:end);
%% grid of log-variance
mean_var=log(v_0)+(r-1/2*eta^2)*Tmat;
std_var=eta^2*Tmat;
a_v=mean_var-2.0*std_var;
b_v=mean_var+1.7*std_var;
n_v=20;                         % total nodes for log-variance.
low_var=exp(a_v);                % boundary for variance
up_var=exp(b_v);

[nodes,w]=lgwt(n_v,a_v,b_v);       % compute corresponding nodes and weights for Guass-Legendre.
Zeta=fliplr(nodes.');  



%% exponential factor
a=1.4;


%% preliminary compution

%% evaluation of Gamma function

gamma1m=2+1i*omega_m;gamma1p=2+1i*omega_p;
gamma2m1=1i*(omega_m-omega_p(1));gamma2m2=1i*(omega_m(1)-omega_p(2:end));
gamma2p1=1i*(omega_p-omega_m(1));gamma2p2=1i*(omega_p(1)-omega_m(2:end));

l=ceil(log2(2*(n_l+n_r+1)-1));
cm=zeros(2^l,1);cp=zeros(2^l,1);

gamma1m1=cgamma(gamma1m(n_l+1:end));
gamma1p1=cgamma(gamma1p(n_l+1:end));
Gamma1m=exp(a*real(omega_m)).*[flipud(conj(gamma1m1(2:n_l+1)));gamma1m1];
Gamma1p=exp(a*real(omega_p)).*[flipud(conj(gamma1p1(2:n_l+1)));gamma1p1];
 
gamma2m11=cgamma(gamma2m1);
gamma2m21=conj(gamma2m11(2:end));
Gamma2m1=exp(a*real(gamma2m1/1i)).*gamma2m11;
Gamma2m2=exp(a*real(gamma2m2/1i)).*gamma2m21;

gamma2p11=cgamma(gamma2p1);
gamma2p21=conj(gamma2p11(2:end));
Gamma2p1=exp(a*real(gamma2p1/1i)).*gamma2p11;
Gamma2p2=exp(a*real(gamma2p2/1i)).*gamma2p21;



Gammainvm=1./Gamma1m(n_l+1:end);
Gammainvp=1./Gamma1p(n_l+1:end);

cm(1:n_l+n_r+1)=Gamma2m1;
cm(2^l-(n_l+n_r-1):2^l)=flipud(Gamma2m2);
cp(1:n_l+n_r+1)=Gamma2p1;
cp(2^l-(n_l+n_r-1):2^l)=flipud(Gamma2p2);




%% product of conditional ChF of increments and transition density of v
%% Hull-White
Chara_m=zeros(n_l+n_r+1,n_v,n_v);
Chara_p=zeros(n_l+n_r+1,n_v,n_v);

Zeta_t=reshape(Zeta,1,[],1);
Zeta_s=reshape(Zeta,1,1,[]);

n_z=300;L_z=4;
delta_z=L_z/n_z;
z_grid=delta_z+(0:delta_z:L_z).';
w_z=[1/2;ones(n_z,1)];
kappa_hat_m=rho*kappa*1i*omegam1./sqrt(-1i*omegam1+(1-rho^2)*omegam1.^2);
kappa_hat_p=rho*kappa*1i*omegap1./sqrt(-1i*omegap1+(1-rho^2)*omegap1.^2);
theta_hat_m=eta*sqrt(-1i*omegam1+(1-rho^2)*omegam1.^2);
theta_hat_p=eta*sqrt(-1i*omegap1+(1-rho^2)*omegap1.^2);

n_xi=25;L_xi=5;
delta_xi=L_xi/n_xi;
xi_grid=(0:delta_xi:L_xi).';
w_xi=[1/2;ones(n_xi,1)];

p_hat_m=0;
p_hat_p=0;
p_hat_0_m=0;
p_hat_0_p=0;
for m=1:n_z+1
    Theta_v_1_m=4*sqrt(exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_m.*z_grid(m))./(eta^2*z_grid(m).*sqrt(pi^3*eta^2*delta_t/2));
    Theta_v_1_p=4*sqrt(exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_p.*z_grid(m))./(eta^2*z_grid(m).*sqrt(pi^3*eta^2*delta_t/2));
    Theta_v0_1_m=4*sqrt(exp(Zeta)/v_0+exp(Zeta).*theta_hat_m.*z_grid(m))./(eta^2*z_grid(m).*sqrt(pi^3*eta^2*delta_t/2));
    Theta_v0_1_p=4*sqrt(exp(Zeta)/v_0+exp(Zeta).*theta_hat_p.*z_grid(m))./(eta^2*z_grid(m).*sqrt(pi^3*eta^2*delta_t/2));
    
    Theta_v_2_m=0;
    Theta_v_2_p=0;
    Theta_v0_2_m=0;
    Theta_v0_2_p=0;
    for n=1:n_xi+1
        Theta_v_2_m=Theta_v_2_m+exp(2*(pi^2-xi_grid(n).^2)/(eta^2*delta_t)-4*sqrt(exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_m.*z_grid(m))./(eta^2*z_grid(m)).*cosh(xi_grid(n))).*...
            sinh(xi_grid(n)).*sin(4*pi*xi_grid(n)/(eta^2*delta_t)).*w_xi(n)*delta_xi;
        Theta_v_2_p=Theta_v_2_p+exp(2*(pi^2-xi_grid(n).^2)/(eta^2*delta_t)-4*sqrt(exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_p.*z_grid(m))./(eta^2*z_grid(m)).*cosh(xi_grid(n))).*...
            sinh(xi_grid(n)).*sin(4*pi*xi_grid(n)/(eta^2*delta_t)).*w_xi(n)*delta_xi;
        Theta_v0_2_m=Theta_v0_2_m+exp(2*(pi^2-xi_grid(n).^2)/(eta^2*delta_t)-4*sqrt(exp(Zeta)/v_0+exp(Zeta).*theta_hat_m.*z_grid(m))./(eta^2*z_grid(m)).*cosh(xi_grid(n))).*...
            sinh(xi_grid(n)).*sin(4*pi*xi_grid(n)/(eta^2*delta_t)).*w_xi(n)*delta_xi;
        Theta_v0_2_p=Theta_v0_2_p+exp(2*(pi^2-xi_grid(n).^2)/(eta^2*delta_t)-4*sqrt(exp(Zeta)/v_0+exp(Zeta).*theta_hat_p.*z_grid(m))./(eta^2*z_grid(m)).*cosh(xi_grid(n))).*...
            sinh(xi_grid(n)).*sin(4*pi*xi_grid(n)/(eta^2*delta_t)).*w_xi(n)*delta_xi;
    end
    p_hat_m=p_hat_m+(exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_m.*z_grid(m)).^((kappa_hat_m-1/2*eta^2)./eta^2).*...
        exp(-(kappa_hat_m-1/2*eta^2).^2/(2*eta^2)*delta_t-2*(1+exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_m.*z_grid(m))./(eta^2*z_grid(m))).*...
        Theta_v_1_m.*Theta_v_2_m./(2*z_grid(m).*exp(Zeta_t)).*exp(Zeta_t).*w_z(m)*delta_z;
    p_hat_p=p_hat_p+(exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_p.*z_grid(m)).^((kappa_hat_p-1/2*eta^2)./eta^2).*...
        exp(-(kappa_hat_p-1/2*eta^2).^2/(2*eta^2)*delta_t-2*(1+exp(Zeta_t-Zeta_s)+exp(Zeta_t).*theta_hat_p.*z_grid(m))./(eta^2*z_grid(m))).*...
        Theta_v_1_p.*Theta_v_2_p./(2*z_grid(m).*exp(Zeta_t)).*exp(Zeta_t).*w_z(m)*delta_z;
    p_hat_0_m=p_hat_0_m+(exp(Zeta)/v_0+exp(Zeta).*theta_hat_m.*z_grid(m)).^((kappa_hat_m-1/2*eta^2)./eta^2).*...
        exp(-(kappa_hat_m-1/2*eta^2).^2/(2*eta^2)*delta_t-2*(1+exp(Zeta)/v_0+exp(Zeta).*theta_hat_m.*z_grid(m))./(eta^2*z_grid(m))).*...
        Theta_v0_1_m.*Theta_v0_2_m./(2*z_grid(m).*exp(Zeta)).*exp(Zeta).*w_z(m)*delta_z;
    p_hat_0_p=p_hat_0_p+(exp(Zeta)/v_0+exp(Zeta).*theta_hat_p.*z_grid(m)).^((kappa_hat_p-1/2*eta^2)./eta^2).*...
        exp(-(kappa_hat_p-1/2*eta^2).^2/(2*eta^2)*delta_t-2*(1+exp(Zeta)/v_0+exp(Zeta).*theta_hat_p.*z_grid(m))./(eta^2*z_grid(m))).*...
        Theta_v0_1_p.*Theta_v0_2_p./(2*z_grid(m).*exp(Zeta)).*exp(Zeta).*w_z(m)*delta_z;
%     m
end

n_1m=-sqrt(-1i*omegam1+(1-rho^2)*omegam1.^2)/eta;
n_1p=-sqrt(-1i*omegap1+(1-rho^2)*omegap1.^2)/eta;
n_2m=-kappa*rho*1i*omegam1./(eta^3*n_1m)-kappa/eta^2;
n_2p=-kappa*rho*1i*omegap1./(eta^3*n_1p)-kappa/eta^2;
phi_vm=exp(n_1m.*(exp(Zeta_s)-exp(Zeta_t))).*exp((Zeta_s-Zeta_t).*n_2m);
phi_vp=exp(n_1p.*(exp(Zeta_s)-exp(Zeta_t))).*exp((Zeta_s-Zeta_t).*n_2p);
C_m=1/2*eta^2*n_2m.^2+(kappa-1/2*eta^2)*n_2m;
C_p=1/2*eta^2*n_2p.^2+(kappa-1/2*eta^2)*n_2p;
Charam1=exp(-1i*omegam1*(r-q)*delta_t-1i*omegam1*rho/eta.*(exp(Zeta_t)-exp(Zeta_s))+C_m*delta_t).*phi_vm.*p_hat_m;
Chara_m(n_l+1:end,:,:)=Charam1;
Chara_m(1:n_l,:,:)=conj(Charam1(n_l+1:-1:2,:,:));
Charap1=exp(-1i*omegap1*(r-q)*delta_t-1i*omegap1*rho/eta.*(exp(Zeta_t)-exp(Zeta_s))+C_p*delta_t).*phi_vp.*p_hat_p;
Chara_p(n_l+1:end,:,:)=Charap1;
Chara_p(1:n_l,:,:)=conj(Charap1(n_l+1:-1:2,:,:));

phi_v0m=exp(n_1m.*(v_0-exp(Zeta))).*(v_0./exp(Zeta)).^n_2m;
phi_v0p=exp(n_1p.*(v_0-exp(Zeta))).*(v_0./exp(Zeta)).^n_2p;
Charavm01=exp(-1i*omegam1*(r-q)*delta_t-1i*omegam1*rho/eta.*(exp(Zeta)-v_0)+C_m*delta_t).*phi_v0m.*p_hat_0_m;
Charavp01=exp(-1i*omegap1*(r-q)*delta_t-1i*omegap1*rho/eta.*(exp(Zeta)-v_0)+C_p*delta_t).*phi_v0p.*p_hat_0_p;
 



%% backward induction

W_hat_m_1=Gamma1m./(omega_m.*(1i-omega_m));
W_hat_p_1=Gamma1p./(omega_p.*(1i-omega_p));

W_hat_m=sum(W_hat_m_1.*Chara_m.*w.',2);
W_hat_p=sum(W_hat_p_1.*Chara_p.*w.',2);
W_hat_m=reshape(W_hat_m,n_l+n_r+1,n_v,1);
W_hat_p=reshape(W_hat_p,n_l+n_r+1,n_v,1);


 
for t_k=1:N_monitor-2
    W_m=W_hat_m;
    W_p=W_hat_p;
    T_fft_m=ifft(fft(cm).*fft([W_p;zeros(2^l-(n_l+n_r+1),n_v)]));
    T_fft_p=ifft(fft(cp).*fft([W_m;zeros(2^l-(n_l+n_r+1),n_v)]));
    T_m=T_fft_m(1:n_l+n_r+1,:);
    T_p=T_fft_p(1:n_l+n_r+1,:);
    W_hat_m=Gamma1m.*(1./(1i*omega_m)-1./(1+1i*omega_m)*(1-exp((r-q)*(t_k+1)*delta_t))/(1-exp((r-q)*delta_t))).*sum(Chara_m.*w.',2)+1/(2*pi)*sum(Chara_m.*T_m.*w.',2)*delta_omega;
    W_hat_p=Gamma1p.*(1./(1i*omega_p)-1./(1+1i*omega_p)*(1-exp((r-q)*(t_k+1)*delta_t))/(1-exp((r-q)*delta_t))).*sum(Chara_p.*w.',2)+sum(Chara_p.*w.'.*(W_p+1/(2*pi)*T_p*delta_omega),2);
    W_hat_m=reshape(W_hat_m,n_l+n_r+1,n_v,1);
    W_hat_p=reshape(W_hat_p,n_l+n_r+1,n_v,1);
end

W_m=W_hat_m;
W_p=W_hat_p;
T_fft_m=ifft(fft(cm).*fft([W_p;zeros(2^l-(n_l+n_r+1),n_v)]));
T_fft_p=ifft(fft(cp).*fft([W_m;zeros(2^l-(n_l+n_r+1),n_v)]));
T_m=T_fft_m(1:n_l+n_r+1,:);
T_p=T_fft_p(1:n_l+n_r+1,:);


x_0=-log((N_monitor+1)*K./S_0-1);

ww=[1/2,ones(1,n_r)];

%% Fourier inversion
W_hat_m_v0=Gamma1m(n_l+1:end).*(1./(1i*omegam1)-1./(1+1i*omegam1)*(1-exp((r-q)*Tmat))/(1-exp((r-q)*delta_t))).*Charavm01*w+1/(2*pi)*(Charavm01.*T_m(n_l+1:end,:))*delta_omega*w;
V_c_m=exp(-r*Tmat)/(N_monitor+1)/pi*real(max((N_monitor+1)*K-S_0,0).*exp(-1i*x_0*omegam1.').*ww*(W_hat_m_v0.*Gammainvm)*delta_omega);

% Re_error_call=abs(V_c_m-V_c_benchmark)./V_c_benchmark;


W_hat_p_v0=Gamma1p(n_l+1:end).*(1./(1i*omegap1)-1./(1+1i*omegap1)*(1-exp((r-q)*Tmat))/(1-exp((r-q)*delta_t))).*Charavp01*w+(Charavp01.*(W_p(n_l+1:end,:)+1/(2*pi)*T_p(n_l+1:end,:)*delta_omega))*w;
V_c_p=exp(-r*Tmat)/(N_monitor+1)/pi*real(max((N_monitor+1)*K-S_0,0).*exp(-1i*x_0*omegap1.').*ww*(W_hat_p_v0.*Gammainvp)*delta_omega);


toc;