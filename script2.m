function script2(deltap,deltam,China)

global sigma gamma b alpha1 alpha2 Mort_Rate pC pN pR wC wR NC scenario scenario_p JDI JAI K Pos Posl N J T1 T0 Country I0 D0 H0 delay Type J_id Z_tau eta;

%% Estimation of alpha
Ka = J-JAI-JDI; A = zeros(Ka,1);
for i = 0:Ka-1
    A(i+1) = (D0(J-i)-D0(J-i-JDI))/(I0(J-i-JDI)-D0(J-i-JDI)-H0(J-i-JDI)-I0(J-i-JDI-JAI+1)+D0(J-i-JDI-JAI+1)+H0(J-i-JDI-JAI+1))/Mort_Rate;
end
alpha = mean(A);
figure('Name',['Alpha for ' Country],'Position',Pos);
plot(1:Ka,A,'b',1:Ka,alpha*ones(Ka,1),'b:');
I0 = alpha1*(I0-D0-H0); R0 = D0+alpha2*H0; % measured recovered 

%% Schedule of quarantine
M = T1+7*(wC+wR)*NC; % number of days to simulate
t = T0 + caldays(0:T1+7*(wC+wR)*NC-1);
r = pC*ones(M,1); r(1:T1) = pN;
r(T1+wC*7+1:T1+(wC+wR)*7) = pR;
ticks = [t(T1+1) t(T1+wC*7) t(T1+(wC+wR)*7)];
for i = 2:NC
    r(T1+(i-1)*(wC+wR)*7+wC*7+1:T1+i*(wC+wR)*7) = pR;
    ticks = [ticks t(T1+(i-1)*(wC+wR)*7+wC*7) t(T1+i*(wC+wR)*7)];
end
figure('Name',['Number of contacts for ' Country],'Position',Pos);
plot(t,r); ylim([pC-1 pN+1]);
ax = gca; set(ax,'xTick',ticks,'XGrid','on');

%% Simulation for Chinese data
if (China == 1)
    b_ = mean(b); gamma_ = mean(gamma);
    R = ones(M,1); R(1) = R0(1);
    I = ones(M,1); I(1) = I0(1);
    E = ones(M,1); E(1) = (I0(2)-(1-gamma_)*I0(1))/sigma; % alpha*I(1);
    S = ones(M,1); S(1) = N-R(1)-I(1)-E(1);

    for i = 1:M-1
        S(i+1) = S(i)-b_*(pC*I(i)+r(delay_ind(i-delay))*E(i))*S(i)/N/(1+eta*I(i));
        E(i+1) = (1-sigma)*E(i)+b_*(pC*I(i)+r(delay_ind(i-delay))*E(i))*S(i)/N/(1+eta*I(i));
        I(i+1) = (1-gamma_)*I(i)+sigma*E(i);
        R(i+1) = R(i)+gamma_*I(i);
    end

    figure('Name',['Simulation with Chinese data for scenario ' num2str(scenario) ' in ' Country],'Position',Pos);
    semilogy(t,E,'r',t,I,'b',t,R,'g',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); legend('E','I','R','Location','southeast');
    ax = gca; set(ax,'xTick',ticks,'XGrid','on');

    figure('Name',['Verification of Chinese parameters by measured data for ' Country],'Position',Pos);
    subplot(2,1,1); semilogy(t(1:J),I(1:J),'b',t(1:J),I0,'b--'); title('I'); xlim([t(1) t(J)]);
    subplot(2,1,2); semilogy(t(1:J),R(1:J),'g',t(1:J),R0,'g--'); title('R'); xlim([t(1) t(J)]);
end

%% Estimation of parameters
E0 = ones(J-J_id-1,1); S0 = ones(J-J_id-1,1);
for i = 2:J-J_id
    E0(i-1) = (I0(i)-I0(i-1)+R0(i)-R0(i-1))/sigma;
%     if (i > 2)&&(E0(i-1) < E0(i-2))
%         E0(i-1) = E0(i-2);
%     elseif (i == 2)&&(E0(1) < 0)
%         E0(1) = 1;
%     end
    S0(i-1) = N-R0(i-1)-I0(i-1)-E0(i-1);
end

% Delay
dI0 = log(I0(2:J-J_id))-log(I0(1:J-J_id-1));
if (J-J_id-2-T1 > 0)
    % Method 1
    E_tau = zeros(J-J_id-2-T1,1);
    for i = T1+1:J-J_id-1
        E_tau(i-T1) = norm(dI0-[mean(dI0(1:i-1))*ones(1,i-1) mean(dI0(i:J-J_id-1))*ones(1,J-J_id-i)]);
    end
    [y, i_min] = min(E_tau); i_min = max(i_min)
    %tau = min(delay,i_min)
    tau = delay;
    % Method 2
    % a1 = (log(E0(Z_tau))-log(E0(1)))/(Z_tau-1); b1 = (Z_tau*log(E0(1))-log(E0(Z_tau)))/(Z_tau-1);
    % a2 = (log(E0(J-J_id-1))-log(E0(J-J_id-Z_tau)))/(Z_tau-1); b2 = ((J-J_id-1)*log(E0(J-J_id-Z_tau))-(J-J_id-Z_tau)*log(E0(J-J_id-1)))/(Z_tau-1);
    % xi = ceil((Z_tau*log(E0(1))-log(E0(Z_tau))-(J-J_id-1)*log(E0(J-J_id-Z_tau))+(J-J_id-Z_tau)*log(E0(J-J_id-1)))/(log(E0(J-J_id-1))-log(E0(J-J_id-Z_tau))-log(E0(Z_tau))+log(E0(1)))-T1)
    ab1 = [1:Z_tau;ones(1,Z_tau)]'\log(E0(1:Z_tau));
    ab2 = [J-J_id-Z_tau:J-J_id-1;ones(1,Z_tau)]'\log(E0(J-J_id-Z_tau:J-J_id-1));
    xi = ceil((ab1(2)-ab2(2))/(ab2(1)-ab1(1)))-T1
else
    tau = delay; ab1 = [1; 1]; ab2 = [-1; 1]; xi = delay;
end
if tau == delay
    i_min = min(J-J_id-1,T1+delay);
end

% gamma & b
if (scenario_p == 1)
    B = zeros(J-J_id-K,1); G = zeros(J-J_id-K,1);
    for k = 0:J-J_id-1-K
        den = 0; num = 0;
        for i = 1:K
            den = den+I0(k+i)^2;
            num = num+(R0(k+i+1)-R0(k+i))*I0(k+i);
        end
        G(k+1) = num/den;
    end
    gamma_i = mean(G)
    for k = 0:J-J_id-1-K
        den = 0; num = 0;
        for i = 2:K
            den = den+(S0(k+i-1)*(pC*I0(k+i-1)+r(delay_ind(k+i-1-tau))*E0(k+i-1))/(1+eta*I0(k+i-1)))^2;
            num = num+(E0(k+i)-(1-sigma)*E0(k+i-1))*S0(k+i-1)*(pC*I0(k+i-1)+r(delay_ind(k+i-1-tau))*E0(k+i-1))/(1+eta*I0(k+i-1));
        end
        B(k+1) = N*num/den;
    end
    b_i = mean(B(B>0))
else
    Ki = J-J_id-K; % number of days starting from which estimate the values of parameters
    B = zeros(Ki+1,1); G = zeros(Ki+1,1);
    for k = Ki:-1:0
        den = 0; num = 0;
        for i = 1:J-J_id-k-1
            den = den+I0(i)^2;
            num = num+(R0(i+1)-R0(i))*I0(i);
        end
        G(Ki-k+1) = num/den;
    end
    gamma_i = mean(G)
    for k = Ki:-1:0
        den = 0; num = 0;
        for i = 2:J-J_id-k-1
            den = den+(S0(i-1)*(pC*I0(i-1)+r(delay_ind(i-1-tau))*E0(i-1))/(1+eta*I0(i-1)))^2;
            num = num+(E0(i)-(1-sigma)*E0(i-1))*S0(i-1)*(pC*I0(i-1)+r(delay_ind(i-1-tau))*E0(i-1))/(1+eta*I0(i-1));
        end
        B(Ki-k+1) = N*num/den;
    end
    b_i = mean(B(B>0))
end
figure('Name',['Identified parameters on a window for ' Country],'Position',Pos);
subplot(2,2,1); plot(t(K+1:J-J_id),G,'m',t(K+1:J-J_id),gamma_i*ones(J-J_id-K,1),'m:'); title('\gamma'); xlim([t(K+1) t(J-J_id)]);
subplot(2,2,2); plot(t(K+1:J-J_id),B,'r',t(K+1:J-J_id),b_i*ones(J-J_id-K,1),'r:'); title('b'); xlim([t(K+1) t(J-J_id)]);
subplot(2,2,3); plot(t(2:J-J_id),dI0,'c',t(2:J-J_id),[mean(dI0(1:i_min-1))*ones(1,i_min-1) mean(dI0(i_min:J-J_id-1))*ones(1,J-J_id-i_min)],'c--'); title('\tau 1 method'); xlim([t(2) t(J-J_id)]);
subplot(2,2,4); plot(t(2:J-J_id),log(E0),'b',t(2:J-J_id),ab1(1)*(1:J-J_id-1)+ab1(2),'b--',t(2:J-J_id),ab2(1)*(1:J-J_id-1)+ab2(2),'b--'); title('\tau 2 method'); xlim([t(2) t(J-J_id)]); ylim([0.9*min(log(E0)) 1.1*max(ab1(1)*(xi+T1)+ab1(2),max(log(E0)))]);

%% Simulation for identified data
R(1) = R0(1); I(1) = I0(1); E(1) = E0(1); S(1) = S0(1);

for i = 1:M-1
    S(i+1) = S(i)-b_i*(pC*I(i)+r(delay_ind(i-tau))*E(i))*S(i)/N/(1+eta*I(i));
    E(i+1) = (1-sigma)*E(i)+b_i*(pC*I(i)+r(delay_ind(i-tau))*E(i))*S(i)/N/(1+eta*I(i));
    I(i+1) = (1-gamma_i)*I(i)+sigma*E(i);
    R(i+1) = R(i)+gamma_i*I(i);
end

figure('Name',['Simulation with identified data for scenario ' num2str(scenario) ' in ' Country],'Position',Pos);
semilogy(t,E,'r',t,I,'b',t,R,'g',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); legend('E','I','R','Location','southeast');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');

figure('Name',['Verification of identified parameters by measured data for ' Country],'Position',Pos);
subplot(2,2,1); semilogy(t(1:J),I(1:J),'b',t(1:J-J_id),I0(1:J-J_id),'b--',t(J-J_id+1:J),I0(J-J_id+1:J),'bs'); title('I'); xlim([t(1) t(J)]);
subplot(2,2,2); semilogy(t(1:J),R(1:J),'g',t(1:J-J_id),R0(1:J-J_id),'g--',t(J-J_id+1:J),R0(J-J_id+1:J),'gs'); title('R'); xlim([t(1) t(J)]);
subplot(2,2,3); semilogy(t(1:J-J_id-1),E(1:J-J_id-1),'r',t(1:J-J_id-1),E0,'r--'); title('E'); xlim([t(1) t(J-J_id-1)]);
subplot(2,2,4); semilogy(t(1:J-J_id-1),S(1:J-J_id-1),'m',t(1:J-J_id-1),S0,'m--'); title('S'); xlim([t(1) t(J-J_id-1)]);
disp(['Model error for I = ' num2str(norm(I(1:J)-I0(1:J)),'%10.3e') ' for R = ' num2str(norm(R(1:J)-R0(1:J)),'%10.3e')]);
disp(' ');

%% Interval observer simulaiton for deviation of all parameters
delta_ = [1-deltam 1+deltap]; % variation for intervals
gamma_ = gamma_i*delta_; b_ = b_i*delta_; r_ = ceil(r.*delta_); sigma_ = sigma*delta_;
Rm = zeros(M,1); RM = zeros(M,1); Im = zeros(M,1); IM = zeros(M,1);
Em = zeros(M,1); EM = zeros(M,1); Sm = zeros(M,1); SM = zeros(M,1);
% inirial conditions
Rm(1) = delta_(1)*R(1); RM(1) = delta_(2)*R(1);
Im(1) = delta_(1)*I(1); IM(1) = delta_(2)*I(1);
Em(1) = delta_(1)*E(1); EM(1) = delta_(2)*E(1);
Sm(1) = delta_(1)*S(1); SM(1) = delta_(2)*S(1);

for i = 1:M-1
    Sm(i+1) = (1-b_(2)*(pC*IM(i)+r_(delay_ind(i-tau),2)*EM(i))/N/(1+eta*Im(i)))*Sm(i);
    Em(i+1) = (1-sigma_(2)+b_(1)*r_(delay_ind(i-tau),1)*Sm(i)/N/(1+eta*IM(i)))*Em(i)+b_(1)*pC*Im(i)*Sm(i)/N/(1+eta*IM(i));
    Im(i+1) = (1-gamma_(2))*Im(i)+sigma_(1)*Em(i);
    Rm(i+1) = Rm(i)+gamma_(1)*Im(i);
    SM(i+1) = min(N,(1-b_(1)*(pC*Im(i)+r_(delay_ind(i-tau),1)*Em(i))/N/(1+eta*IM(i)))*SM(i));
    EM(i+1) = min(N,(1-sigma_(1)+b_(2)*r_(delay_ind(i-tau),2)*SM(i)/N/(1+eta*Im(i)))*EM(i)+b_(2)*pC*IM(i)*SM(i)/N/(1+eta*Im(i)));
    IM(i+1) = min(N,(1-gamma_(1))*IM(i)+sigma_(2)*EM(i));
    RM(i+1) = min(N,RM(i)+gamma_(2)*IM(i));
end

figure('Name',['Simulation of interval predictor for deviation of all parameters for ' Country],'Position',Posl);
subplot(2,2,1); semilogy(t,E,'r',t,Em,'r:',t,EM,'r--',t(1:J-J_id-1),E0,'ro',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('E');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');
subplot(2,2,2); semilogy(t,I,'b',t,Im,'b:',t,IM,'b--',t(1:J-J_id),I0(1:J-J_id),'bo',t(J-J_id+1:J),I0(J-J_id+1:J),'bs',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('I');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');
subplot(2,2,3); semilogy(t,R,'g',t,Rm,'g:',t,RM,'g--',t(1:J-J_id),R0(1:J-J_id),'go',t(J-J_id+1:J),R0(J-J_id+1:J),'gs',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('R');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');
subplot(2,2,4); semilogy(t,S,'m',t,Sm,'m:',t,SM,'m--',t(1:J-J_id-1),S0,'mo',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('S');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');

figure('Name',['Interval predictor of infectives for deviation of all parameters for ' Country],'Position',Pos);
if (Type == 1)
    semilogy(t,Im,'b:',t,IM,'b--',t(1:J-J_id),I0(1:J-J_id),'bo',t(J-J_id+1:J),I0(J-J_id+1:J),'bs',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('I'); ylim([min(Im) max(IM)]);
else
    plot(t(1:J-J_id),I0(1:J-J_id),'bo',t(J-J_id+1:J),I0(J-J_id+1:J),'bs'); xlim([t(1) t(M)]); title('I'); %,t,N*ones(M,1),'k:'
end
ax = gca; set(ax,'xTick',ticks,'XGrid','on'); hold on

% initial conditions from the last date used for identification
Rm(J-J_id-1) = delta_(1)*R0(J-J_id-1); RM(J-J_id-1) = delta_(2)*R0(J-J_id-1);
Im(J-J_id-1) = delta_(1)*I0(J-J_id-1); IM(J-J_id-1) = delta_(2)*I0(J-J_id-1);
Em(J-J_id-1) = delta_(1)*E0(J-J_id-1); EM(J-J_id-1) = delta_(2)*E0(J-J_id-1);
Sm(J-J_id-1) = delta_(1)*S0(J-J_id-1); SM(J-J_id-1) = delta_(2)*S0(J-J_id-1);
R_ = zeros(M,1); E_ = zeros(M,1); I_ = zeros(M,1); S_ = zeros(M,1);
E_(J-J_id-1) = E0(J-J_id-1); R_(J-J_id-1) = R0(J-J_id-1); I_(J-J_id-1) = I0(J-J_id-1); S_(J-J_id-1) = S0(J-J_id-1);

for i = J-J_id-1:M-1
    S_(i+1) = S_(i)-b_i*(pC*I_(i)+r(delay_ind(i-tau))*E_(i))*S_(i)/N/(1+eta*I(i));
    E_(i+1) = (1-sigma)*E_(i)+b_i*(pC*I_(i)+r(delay_ind(i-tau))*E_(i))*S_(i)/N/(1+eta*I(i));
    I_(i+1) = (1-gamma_i)*I_(i)+sigma*E_(i);
    R_(i+1) = R_(i)+gamma_i*I_(i);
    Sm(i+1) = (1-b_(2)*(pC*IM(i)+r_(delay_ind(i-tau),2)*EM(i))/N/(1+eta*Im(i)))*Sm(i);
    Em(i+1) = (1-sigma_(2)+b_(1)*r_(delay_ind(i-tau),1)*Sm(i)/N/(1+eta*IM(i)))*Em(i)+b_(1)*pC*Im(i)*Sm(i)/N/(1+eta*IM(i));
    Im(i+1) = (1-gamma_(2))*Im(i)+sigma_(1)*Em(i);
    Rm(i+1) = Rm(i)+gamma_(1)*Im(i);
    SM(i+1) = min(N,(1-b_(1)*(pC*Im(i)+r_(delay_ind(i-tau),1)*Em(i))/N/(1+eta*IM(i)))*SM(i));
    EM(i+1) = min(N,(1-sigma_(1)+b_(2)*r_(delay_ind(i-tau),2)*SM(i)/N/(1+eta*Im(i)))*EM(i)+b_(2)*pC*IM(i)*SM(i)/N/(1+eta*Im(i)));
    IM(i+1) = min(N,(1-gamma_(1))*IM(i)+sigma_(2)*EM(i));
    RM(i+1) = min(N,RM(i)+gamma_(2)*IM(i));
end

if (Type == 1)
    semilogy(t(J-J_id-1:M),I_(J-J_id-1:M),'r',t(J-J_id-1:M),Im(J-J_id-1:M),'b:',t(J-J_id-1:M),IM(J-J_id-1:M),'b--','LineWidth',2);
else
    plot(t(J-J_id-1:M),I_(J-J_id-1:M),'r',t(J-J_id-1:M),Im(J-J_id-1:M),'b:','LineWidth',2); %,t(J-J_id-1:M),IM(J-J_id-1:M),'b--'
end
[y,x] = max(I_(J-J_id-1:M)); 
if (Type == 1)
    s = 1.6*y;
else
    s = 1.05*y;
end
if (x < M-J-J_id+2)
    x = t(x+J-J_id-2);
    line('XData', [x x t(1)],'YData',[1 y y],'LineWidth',1,'LineStyle', ':','Color','r');
    text(x, s, ['I_{peak} = ' num2str(floor(y)) ' on ' datestr(x,'dd-mmm')]);
end
[y,x] = max(Im(J-J_id-1:M));
if (Type == 1)
    s = 1.7*y;
else
    s = 1.1*y;
end
if (x < M-J-J_id+2)
    x = t(x+J-J_id-2);
    line('XData', [x x t(1)],'YData',[1 y y],'LineWidth',1,'LineStyle', ':','Color','r');
    text(x, s, ['I_{peak}^{optimiste} = ' num2str(floor(y)) ' on ' datestr(x,'dd-mmm')]);
end

%% Interval observer simulation for deviation of sigma
gamma_ = gamma_i*[1 1]; b_ = b_i*[1 1]; r_ = ceil(r.*[1 1]); sigma_ = [1/12 1/2];
Rm = zeros(M,1); RM = zeros(M,1); Im = zeros(M,1); IM = zeros(M,1);
Em = zeros(M,1); EM = zeros(M,1); Sm = zeros(M,1); SM = zeros(M,1);
% inirial conditions
Rm(1) = delta_(1)*R(1); RM(1) = delta_(2)*R(1);
Im(1) = delta_(1)*I(1); IM(1) = delta_(2)*I(1);
Em(1) = delta_(1)*E(1); EM(1) = delta_(2)*E(1);
Sm(1) = delta_(1)*S(1); SM(1) = delta_(2)*S(1);

for i = 1:M-1
    Sm(i+1) = (1-b_(2)*(pC*IM(i)+r_(delay_ind(i-tau),2)*EM(i))/N/(1+eta*Im(i)))*Sm(i);
    Em(i+1) = (1-sigma_(2)+b_(1)*r_(delay_ind(i-tau),1)*Sm(i)/N/(1+eta*IM(i)))*Em(i)+b_(1)*pC*Im(i)*Sm(i)/N/(1+eta*IM(i));
    Im(i+1) = (1-gamma_(2))*Im(i)+sigma_(1)*Em(i);
    Rm(i+1) = Rm(i)+gamma_(1)*Im(i);
    SM(i+1) = min(N,(1-b_(1)*(pC*Im(i)+r_(delay_ind(i-tau),1)*Em(i))/N/(1+eta*IM(i)))*SM(i));
    EM(i+1) = min(N,(1-sigma_(1)+b_(2)*r_(delay_ind(i-tau),2)*SM(i)/N/(1+eta*Im(i)))*EM(i)+b_(2)*pC*IM(i)*SM(i)/N/(1+eta*Im(i)));
    IM(i+1) = min(N,(1-gamma_(1))*IM(i)+sigma_(2)*EM(i));
    RM(i+1) = min(N,RM(i)+gamma_(2)*IM(i));
end

figure('Name',['Simulation of interval predictor for deviation of sigma for ' Country],'Position',Posl);
subplot(2,2,1); semilogy(t,E,'r',t,Em,'r:',t,EM,'r--',t(1:J-J_id-1),E0,'ro',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('E');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');
subplot(2,2,2); semilogy(t,I,'b',t,Im,'b:',t,IM,'b--',t(1:J-J_id),I0(1:J-J_id),'bo',t(J-J_id+1:J),I0(J-J_id+1:J),'bs',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('I');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');
subplot(2,2,3); semilogy(t,R,'g',t,Rm,'g:',t,RM,'g--',t(1:J-J_id),R0(1:J-J_id),'go',t(J-J_id+1:J),R0(J-J_id+1:J),'gs',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('R');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');
subplot(2,2,4); semilogy(t,S,'m',t,Sm,'m:',t,SM,'m--',t(1:J-J_id-1),S0,'mo',t,N*ones(M,1),'k:'); xlim([t(1) t(M)]); title('S');
ax = gca; set(ax,'xTick',ticks,'XGrid','on');