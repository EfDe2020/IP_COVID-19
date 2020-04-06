close all; clearvars;

global sigma gamma b alpha1 alpha2 Mort_Rate pC pN pR wC wR NC scenario scenario_p JDI JAI K Pos Posl N J T1 T0 Country I0 D0 H0 delay Type J_id;

%% Initialization
% Model parameters
sigma = 1/7; % incubation rate [1/12 1/2]
gamma = [0.0721 0.238]; % rate of recovery or death (for China)
b = [0.05068 0.05429]; % rate of trasmsmission during a contact (for China)
alpha1 = 5; alpha2 = alpha1; % ration between measured/unmeasured and/or between I and E, it is in [5 20]
Mort_Rate = 0.034; % mortality rate is 3.4 (for China)
delay = 10; % nominal value

% Confinement parameters
pC = 2; % number of contacts in quarantin (3 for China)
pN = 12; % number of contacts in normal state (15 for China)
pR = 6; % number of contacts in relaxed quarantin (10 for China)

% Scenario of quarantine 
scenario = 1;
if (scenario == 1)
    wC = 6; % weeks in quarantine
    wR = 2; % weeks in relaxed quarantine
    NC = 3; % number cycles to simulate
else
    wC = 12; % weeks in quarantine
    wR = 4; % weeks in relaxed quarantine
    NC = 2; % number cycles to simulate
end

% Identification and prediciton
delta = .1; % the level of uncertainty for IP
JDI = 5; JAI = 12; % days used for identification of alpha
scenario_p = 1; % the algorithm for identification
if (scenario_p == 1)
    K = 6; % length of moving window to estimate the values of parameters
else
    K = 10; % minimal length of window to estimate the values of parameters
end
J_id = 0; % the number of days used for the model verification

% Misc
Pos = [200 200 900 600];
Posl = [200 200 1500 600];
Type = 1; % type of the plot
Imperial = 0; % 1 = use of Imperial College parameters; 2 = use of Cmmdi parameters
China = 0; % validate chinese parameters or not

%% France
Country = 'France';

% ARS statistics from 12/03
I0 = [2876 3661 4499 5423 6633 7730 9134 10995 12612 14459 16018 19856 22302 25230 29155 32964 37575 40174 44550 52128 56989 59105 64338 68605 70478]; % detected infected
D0 = [61 79 91 127 148 175 244 372 450 562 674 860 1100 1331 1696 1995 2314 2606 3024 3523 4032 4503 5091 5532 5889]; % detected deaths
H0 = [12 12 12 12 12 12 12 12 12 12 2200 2200 3288 3900 4948 5700 6624 7132 7923 9444 10935 12428 14008 15438 16183]; % detected recovered
N = 67064000; % population in FR
J = length(I0);
T1 = 4; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 53; alpha2 = alpha1;
    % Confinement parameters for France
    pC = 1.5; % number of contacts in quarantine
    pN = 4; % number of contacts in normal state 
    pR = 2.75; % number of contacts in relaxed quarantine
elseif (Imperial == 2) % Cmmdi
     alpha1 = 15; alpha2 = alpha1;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script1(delta,0.5*delta,China);
input('Press ENTER to proceed, please')
close all;

%% Italy
Country = 'Italy';

% Statistics from 05/03 (eficiens + hopkins)
I0 = [3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59158 63927 69176 74386 80589 86498 92472 97689 101739 105792 110574 115242 119827 124632 128948]; % detected infected
D0 = [148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215 9134 10023 10779 11591 12428 13155 13915 14681 15362 15887]; % detected deaths
H0 = [414 523 589 622 724 1004 1045 1258 1439 1439 1439 2749 2941 4025 4440 5129 6072 7024 7432 8326 9362 10361 10950 12384 13030 14620 15729 16847 18278 19388 20996 21815]; % detected recovered
N = 60359546; % population in IT
J = length(I0);
T1 = 2; % the date of quarantine
T0 = datetime(2020,3,05); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 64; alpha2 = alpha1;
    % Confinement parameters for Italy
    pC = 1.5; % number of contacts in quarantine
    pN = 3.8; % number of contacts in normal state 
    pR = 2.7; % number of contacts in relaxed quarantine
elseif (Imperial == 2) % Cmmdi
    alpha1 = 17; alpha2 = alpha1;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script1(delta,delta,China);
input('Press ENTER to proceed, please')
close all;

%% Spain
Country = 'Spain';

% Statistics from 12/03 (https://es.statista.com/estadisticas/1104279/personas-con-covid-19-recuperadas-por-dia-espana/ + Hopkins)
I0 = [2950 4209 5753 7753 9383 11273 13716 17147 19980 24926 28572 33089 39675 47600 56188 64059 73235 80110 87956 95923 104118 112065 119199 126168 131646]; % detected infected
D0 = [84 118 136 288 309 497 598 767 982 1326 1720 2182 2696 3434 4365 5138 5982 6803 7716 8464 9387 10348 11198 11947 12641]; % detected deaths
H0 = [189 193 517 517 530 1028 1081 1107 1588 2125 2575 3355 3794 5367 7015 9357 12285 14709 16780 19259 22647 26743 30513 34219 38080]; % detected recovered
N = 46600396; % population in ES
J = length(I0);
T1 = 2; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 95; alpha2 = alpha1;
    % Confinement parameters for Spain
    pC = 1.7; % number of contacts in quarantine
    pN = 4.8; % number of contacts in normal state 
    pR = 3.2; % number of contacts in relaxed quarantine
elseif (Imperial == 2) % Cmmdi
     alpha1 = 18.5; alpha2 = alpha1;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script1(0.5*delta,delta,China);
input('Press ENTER to proceed, please')
close all;

%% Germany
Country = 'Germany';

% Statistics from 12/03 (Hopkins)
I0 = [2369 3062 3795 4838 6012 7156 8198 10999 13957 16662 18610 22672 27436 31554 36508 42288 48582 52547 57298 61913 67366 77981 84794 91159 96092 100123]; % detected infected
D0 = [5 5 8 12 12 12 12 20 31 47 55 86 114 149 198 253 325 389 455 583 775 931 1107 1275 1444 1584]; % detected deaths
H0 = [25 25 46 46 46 67 67 105 115 180 209 266 453 3290 3547 5673 6658 8481 9211 13500 16100 18700 22440 24575 26400 28700]; % detected recovered
N = 83042200; % population in GE (https://en.wikipedia.org/wiki/Demographics_of_Germany)
J = length(I0);
T1 = 10; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 8; alpha2 = alpha1;
    % Confinement parameters for Germany
    pC = 1.8; % number of contacts in quarantine
    pN = 4.6; % number of contacts in normal state 
    pR = 3.2; % number of contacts in relaxed quarantine
elseif (Imperial == 2) % Cmmdi
     alpha1 = 2.13; alpha2 = alpha1;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script1(2*delta,0.5*delta,China);
input('Press ENTER to proceed, please')
close all;

%% Brazil
Country = 'Brazil';

% Statistics from 12/03 (hopkins)
I0 = [77 151 151 200 234 346 529 640 970 1178 1546 1924 2247 2554 2985 3477 3904 4256 4661 5812 6931 8066 9216 10360 11281]; % detected infected
D0 = [0 0 0 0 0 1 4 7 11 18 25 34 46 59 77 93 114 136 165 202 244 327 365 445 487]; % detected deaths
H0 = [0 1 1 1 2 2 2 2 2 2 2 2 2 2 6 6 6 6 127 127 127 127 127 127 127]; % detected recovered
N = 212559417; % population in BR
J = length(I0);
T1 = 2; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % No Imperial
    alpha1 = 5; alpha2 = alpha1;
    % Confinement parameters for Brazil
    pC = 2; % number of contacts in quarantine
    pN = 12; % number of contacts in normal state 
    pR = 6; % number of contacts in relaxed quarantine
elseif (Imperial == 2) % Cmmdi
     alpha1 = 9; alpha2 = alpha1;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script1(delta,delta,China);
input('Press ENTER to proceed, please')
close all;

%% Russia
Country = 'Russia';

% Statistics from 12/03 (Wikipedia + hopkins)
I0 = [45 59 63 93 114 147 199 253 306 367 438 495 658 840 1036 1264 1534 1836 2337 2777 3548 4149 4731 4874 5410]; % detected infected
D0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 2 3 4 8 9 17 24 30 34 43 45 47]; % detected deaths
H0 = [3 3 3 4 5 5 5 12 16 16 17 22 29 38 45 49 64 66 121 190 235 281 333 355 395]; % detected recovered
N = 146745098; % population in RU
J = length(I0);
T1 = 15; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % No Imperial
    alpha1 = 5; alpha2 = alpha1;
    % Confinement parameters for Russia
    pC = 2; % number of contacts in quarantine
    pN = 12; % number of contacts in normal state 
    pR = 6; % number of contacts in relaxed quarantine
elseif (Imperial == 2) % Cmmdi
    alpha1 = 2.5; alpha2 = alpha1;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script1(delta,delta,China);
disp('Stop execution')