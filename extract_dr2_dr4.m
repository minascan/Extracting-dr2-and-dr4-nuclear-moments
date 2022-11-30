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
%F1(1)=-1.849792025703267*(10^5); F2(1)=2.425446341436609*(10^2); F3(1)=-0.635925306138426; F4(1)=0.001037472142435; %LITHIUM-like
% Transition 10
%F1(2)=-0.072740620581839*(10^5); F2(2)=0.090085521902942*(10^2); F3(2)=-0.023664207999006; F4(2)=0.000038259925557; %LITHIUM-like

F1(1)=-1.573512465813949*(10^5); F2(1)= 2.075124708925421*(10^2); F3(1)=-0.545067800809402; F4(1)= 0.000896047878353;
F1(2)= 2.260087431698599*(10^5); F2(2)=-2.939843311017920*(10^2); F3(2)= 0.772387440181695; F4(2)=-0.001266437615842;

% The dr2 and dr4 values used to produce the pseudo-experimental data or
% else the "exact" values
dr2_exp=-0.1638;
dr4_exp=-13.7693;

%The psuedo-experimental data
%nu=nu_exp(:);
%nu=[27422.148184519512; 1084.9898508226213;]; %LITHIUM-like
nu=[23407.79512057857; -33676.28639191137;]; %BERYLLIUM-like


% Here its assumed that the experimental errors are in the 4th decimal
er=zeros(k_max);

for k=1:k_max
er(k) = nu(k) * 10^(-3);  % error for transitions
end

sigma_x = zeros(k_max);

for k=1:k_max
%nu_error(k,:) = [-er(k),0,er(k)];
sigma_x(k,k)  = er(k)^2;
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
K
r  = mldivide(K,nu);

% different error estimate
%sigma_f = sqrt( inv(K) * sigma_x * transpose(inv(K)));

% different error estimate
% pseudoinverse %
Kp = inv(transpose(K)*K) * transpose(K);
sigma_f = sqrt( Kp * sigma_x * transpose(Kp));


T =sprintf('original sum: <dr^2> = %5.4f (%5.4f), <dr^4>= %5.4f (%5.4f)',...
    r(1),sigma_f(1,1),r(2),sigma_f(2,2));
disp(T)
disp(' ')

% -------- Using the y-functions -------------

K1=zeros(k_max,2);

% K1 * y = nu


K1(:,1) = 0.288554*A^(2/3)*F1(:) + 0.350673*A^(4/3)*F2(:)... 
        + 0.448303*A^2*F3(:) + 0.592709*A^(8/3)*F4(:) ;
    
K1(:,2) = 0.0799258*A^(4/3)*F2(:) + 0.172916*A^2*F3(:)...
        + 0.2972*A^(8/3)*F4(:);


y = mldivide(K1,nu);


% Convert back to r2 and r4

% K2 * r = y

K2=zeros(2,2);

K2(1,:) = [3.46556/A^(2/3), 0];
K2(2,:) = [-15.2051/A^(2/3), 12.5116/A^(4/3)];

r  = mldivide(K2,y);

% different error estimate
% pseudoinverse %
Kp = inv(transpose(K1*K2)*(K1*K2)) * transpose(K1*K2);
sigma_f = sqrt( Kp * sigma_x * transpose(Kp));

%sigma_f = sqrt( inv(K2)*inv(K1) * sigma_x * transpose(inv(K2)*inv(K1)));

T =sprintf('y-method sum: <dr^2> = %5.4f (%5.4f), <dr^4>= %5.4f (%5.4f)',...
    r(1),sigma_f(1,1),r(2),sigma_f(2,2));

disp(T)
disp(' ')

% -------- Correct result -------------

r=[dr2_exp;dr4_exp];
disp(' ')
disp('------------------------------------------------------------------- ')
T=sprintf('exact res   : <dr^2> = %5.4f,          <dr^4>= %5.4f',r(1),r(2));
disp(T)
disp('------------------------------------------------------------------- ')
disp(' ')