clear
format long
% average massnumber
global A af cf A1 A2

A2=236;
A1=238;
A=(A2+A1)/2;

% number of transitions that we have
k_max=2;

% Transition 1
F1(1)=-1.849792025703267*(10^5); F2(1)=2.425446341436609*(10^2); F3(1)=-0.635925306138426; F4(1)=0.001037472142435;
% Transition 10
F1(2)=-0.072740620581839*(10^5); F2(2)=0.090085521902942*(10^2); F3(2)=-0.023664207999006; F4(2)=0.000038259925557;

% The dr2 and dr4 values used to produce the pseudo-experimental data or
% else the "exact" values
dr2_exp=-0.1638;
dr4_exp=-13.7693;

%The psuedo-experimental data
%nu=nu_exp(:);
nu=[27422.148184519512; 1084.9898508226213;];

% Here its assumed that the experimental errors are in the 4th decimal
er=zeros(k_max);

for k=1:k_max
er(k) = nu(k) * 10^(-3);  % error for transitions
end

sigma_x = zeros(k_max);

for k=1:k_max
%nu_error(k,:) = [-er(k),0,er(k)];
sigma_x(k,k)  = er(k)^2
end

%sigma_x = [er(1)^2,0; 0, er(2)^2];

disp(' ')
disp('------------------------------------------------------------------- ')
disp('---<dr2> and <dr4>------------------------------------------------- ')
disp(' ')
% -------- Using the r-functions -------------
K=zeros(k_max,2);

% K * r = nu
K(:,1) = F1(:);
K(:,2) = F2(:);

%cond(K,2)



r  = mldivide(K,nu)

K

% different error estimate
INV = inv(K)
III = K * INV

CPX = INV * sigma_x
sigma_f = sqrt( inv(K) * sigma_x * transpose(inv(K)));

% different error estimate
% pseudoinverse %
%KTK_INV = inv(transpose(K)*K)
%Kp = inv(transpose(K)*K) * transpose(K)
%CPX = Kp * sigma_x
%sigma_f = sqrt( Kp * sigma_x * transpose(Kp));


T =sprintf('original sum: <dr^2> = %5.4f (%5.4f), <dr^4>= %5.4f (%5.4f)',...
    r(1),sigma_f(1,1),r(2),sigma_f(2,2));
disp(T)
disp(' ')