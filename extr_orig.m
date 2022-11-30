clear

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
nu=[27422.148184519512; 1084.9898508226213;]


% Here its assumed that the experimental errors are in the 4th decimal

er=zeros(k_max);

for k=1:k_max
er(k) = nu(k) * 10^(-3);  % error for transitions
end

sigma_x = zeros(k_max);

for k=1:k_max
nu_error(k,:) = [-er(k),0,er(k)];
sigma_x(k,k)  = er(k)^2;
end


%sigma_x = [er(1)^2,0; 0, er(2)^2];

e_vec = eye(k_max,k_max);
nu_e = zeros(k_max,k_max^3);

k_e=0;
for k1=1:3
  for k2=1:3
   k_e=k_e + 1;

    nu_e(:,k_e) = nu + e_vec(:,1)*nu_error(1,k1) + e_vec(:,2)*nu_error(2,k2) ;

  end
end
k_e_max = k_e;

disp(' ')
disp('------------------------------------------------------------------- ')
disp('---<dr2> and <dr4>------------------------------------------------- ')
disp(' ')
% -------- Using the r-functions -------------
K=zeros(k_max,2);

% K * r = nu

K(:,1) = F1(:);
K(:,2) = F2(:);

r  = mldivide(K,nu);

% checking the error
error=zeros(2);
for k_e=1:k_e_max
 r_e = mldivide(K,nu_e(:,k_e));
for k=1:2
  error(k) = max(error(k),[abs(r(k)-r_e(k))]);
end
end


% different error estimate
%sigma_f = sqrt( inv(K) * sigma_x * transpose(inv(K)));

% different error estimate
% pseudoinverse %
Kp = inv(transpose(K)*K) * transpose(K);
sigma_f = sqrt( Kp * sigma_x * transpose(Kp));


T =sprintf('original sum: <dr^2> = %5.4f (%5.4f,%5.4f), <dr^4>= %5.4f (%5.4f,%5.4f)',...
    r(1),error(1),sigma_f(1,1),r(2),error(2),sigma_f(2,2));
disp(T)
disp(' ')
