close all; clearvars;

global sigma alpha1 alpha2 Mort_Rate pQ pC pN pR wC wR NC scenario scenario_p JDI JAI K Pos Posl N J T1 T0 Country I0 D0 H0 delay delay_pc Type J_id Z_tau JZ1 JZ2 kappa;

%% Initialization
% Model parameters
sigma = 1/7; % incubation rate [1/12 1/2]
Mort_Rate = 0.034; % mortality rate is 3.4 (for China)
delay = 10; % nominal value
delay_pc = 12; % delay before increasing the stricteness of confinement to pC
kappa = 0.01; % ration between sigma and recovery rate for E
script = @script1;

% Confinement parameters
pQ = 2; % number of contacts in quarantine (3 for China)
pN = 12; % number of contacts in normal state (15 for China)
pR = 6; % number of contacts in relaxed quarantine (10 for China)
pC = 0; % number of contacts in clastering

% Scenario of quarantine 
scenario = 1;
if (scenario == 1)
    wC = 8; % weeks in quarantine
    wR = 4; % weeks in relaxed quarantine
    NC = 3; % number cycles to simulate
else
    wC = 12; % weeks in quarantine
    wR = 4; % weeks in relaxed quarantine
    NC = 2; % number cycles to simulate
end

% Identification and prediction
delta = .1; % the level of uncertainty for IP
JDI = 5; JAI = 12; % days used for identification of alpha
Z_tau = 10; % days for identification of delay
JZ1 = 25; % number of days to use for identification of delay
JZ2 = 44; % number of days to use for identification of delay_pc
scenario_p = 1; % the algorithm for identification
if (scenario_p == 1)
    K = 5; % length of moving window to estimate the values of parameters
else
    K = 10; % minimal length of window to estimate the values of parameters
end
J_id = 0; % the number of days used for the model verification

% Misc
Pos = [200 200 900 600];
Posl = [200 200 1500 600];
Type = 1; % type of the plot
Imperial = 0; % 1 = use of Imperial College parameters; if not use of Cmmdi parameters + delay tuning for each country

%% France
Country = 'France';

% ARS statistics from 12/03
I0 = [2876 3661 4499 5423 6633 7730 9134 10995 12612 14459 16018 19856 22302 25230 29155 32964 37575 40174 44550 52128 56989 59105 64338 68605 70478 74389 78167 82048 86334 90676 93790 95403 98076 103573 106206 108847 109252 111821 112606 114657 117324 119151 120804 122577 124114 124575 128339]; % detected infected
D0 = [63 79 91 127 148 175 264 372 450 562 674 860 1100 1331 1696 1995 2314 2606 3024 3523 4403 5387 6507 7560 8078 8911 10328 10869 12210 13197 13832 14393 14967 15729 17167 17920 18681 19323 19718 20265 20796 21340 21856 22245 22614 22856 23293]; % detected deaths
H0 = [12 12 12 12 12 12 12 12 12 12 2200 2200 3288 3900 4948 5700 6624 7132 7923 9444 10935 12428 14008 15438 16183 17250 19337 21254 23200 24932 26391 27186 27718 28805 30995 32812 34420 35983 36578 37409 39181 40657 42088 43493 44594 44903 45513]; % detected recovered

N = 67064000; % population in FR
J = length(I0);
T1 = 4; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 53; alpha2 = alpha1;
    pQ = 1.5; pN = 4; pR = 2.75;
else
    alpha1 = 23; alpha2 = alpha1;
    delay = 10; delay_pc = 20;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(delta,delta);
input('Press ENTER to proceed, please')
close all;

%% Italy
Country = 'Italy';

% Statistics from 05/03 (eficiens + hopkins)
I0 = [3858 4636 5883 7375 9172 10149 12462 15113 17660 21157 24747 27980 31506 35713 41035 47021 53578 59158 63927 69176 74386 80589 86498 92472 97689 101739 105792 110574 115242 119827 124632 128948 132547 135586 139422 143626 147577 152271 156363 159516 162488 165155 168941 172434 175925 178972 181228 183957 187327 189973 192994 195351 197675 199414]; % detected infected
D0 = [148 197 233 366 463 631 827 1016 1266 1441 1809 2158 2503 2978 3405 4032 4825 5476 6077 6820 7503 8215 9134 10023 10779 11591 12428 13155 13915 14681 15362 15887 16523 17127 17669 18279 18849 19468 19899 20465 21067 21645 22170 22745 23227 23660 24114 24648 25085 25549 25969 26384 26644 26977]; % detected deaths
H0 = [414 523 589 622 724 1004 1045 1258 1439 1439 1439 2749 2941 4025 4440 5129 6072 7024 7432 8326 9362 10361 10950 12384 13030 14620 15729 16847 18278 19388 20996 21815 22837 24392 26491 28470 30455 32534 34211 35435 37130 38092 40164 42727 44927 47055 48877 51600 54543 57576 60498 63120 64928 66624]; % detected recovered

N = 60359546; % population in IT
J = length(I0);
T1 = 2; % the date of quarantine
T0 = datetime(2020,3,05); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 64; alpha2 = alpha1;
    pQ = 1.5; pN = 3.8; pR = 2.7;
else
    alpha1 = 16.7; alpha2 = alpha1;
    delay = 7; delay_pc = 30;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(delta,delta);
input('Press ENTER to proceed, please')
close all;

%% Spain
Country = 'Spain';

% Statistics from 12/03 (https://es.statista.com/estadisticas/1104279/personas-con-covid-19-recuperadas-por-dia-espana/ + Hopkins)
I0 = [2950 4209 5753 7753 9383 11273 13716 17147 19980 24926 28572 33089 39675 47600 56188 64059 73235 80110 87956 95923 104118 112065 119199 126168 131646 136675 141942 148220 153222 158273 166019 166831 170099 174060 180659 184948 190839 194416 198674 200210 204178 208389 213024 219764 223759 226629 229422]; % detected infected
D0 = [84 118 136 288 309 497 598 767 982 1326 1720 2182 2696 3434 4365 5138 5982 6803 7716 8464 9387 10348 11198 11947 12641 13341 14045 14792 15447 16081 16972 17209 17756 18255 18812 19315 20002 20639 20453 20852 21282 21717 22157 22524 22902 23190 23521]; % detected deaths
H0 = [189 193 517 517 530 1028 1081 1107 1588 2125 2575 3355 3794 5367 7015 9357 12285 14709 16780 19259 22647 26743 30513 34219 38080 40437 43208 48021 52165 55668 62391 62391 64727 67504 70853 74797 74797 74797 77357 80587 82514 85915 89250 92355 95708 117727 120832]; % detected recovered

N = 46600396; % population in ES
J = length(I0);
T1 = 2; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 95; alpha2 = alpha1;
    pQ = 1.7; pN = 4.8; pR = 3.2;
else
    alpha1 = 15.3; alpha2 = alpha1;
    delay = 7; delay_pc = 18;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(1.5*delta,0.5*delta);
input('Press ENTER to proceed, please')
close all;

%% Germany
Country = 'Germany';

% Statistics from 12/03 (Hopkins)
I0 = [2369 3062 3795 4838 6012 7156 8198 10999 13957 16662 18610 22672 27436 31554 36508 42288 48582 52547 57298 61913 67366 77981 84794 91159 96092 100123 103375 107663 113296 118235 122171 125452 127854 130072 132210 134753 137698 141397 143724 145742 147065 148453 150648 153129 155054 156513 157770 158758]; % detected infected
D0 = [5 5 8 12 12 12 12 20 31 47 55 86 114 149 198 253 325 389 455 583 775 931 1107 1275 1444 1584 1810 2016 2349 2607 2736 2871 3022 3194 3495 3804 4052 4352 4538 4642 4862 5086 5315 5575 5767 5877 5976 6126]; % detected deaths
H0 = [25 25 46 46 46 67 67 105 115 180 209 266 453 3290 3547 5673 6658 8481 9211 13500 16100 18700 22440 24575 26400 28700 36081 36081 46300 52407 53913 57400 64300 68200 72600 77000 81800 83114 88000 91500 95200 99400 103300 106800 109800 109800 114500 114500]; % detected recovered

N = 83042200; % population in GE (https://en.wikipedia.org/wiki/Demographics_of_Germany)
J = length(I0);
T1 = 10; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % Imperial
    alpha1 = 8; alpha2 = alpha1;
    pQ = 1.8; pN = 4.6; pR = 3.2;
else
    alpha1 = 3.5; alpha2 = alpha1;
    delay = 3; delay_pc = 12;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(delta,delta);
input('Press ENTER to proceed, please')
close all;

%% Brazil
Country = 'Brazil';

% Statistics from 12/03 (hopkins)
I0 = [77 151 151 200 234 346 529 640 970 1178 1546 1924 2247 2554 2985 3477 3904 4256 4661 5812 6931 8066 9216 10360 11281 12232 14049 16188 18176 19943 20964 22318 23723 25684 28912 30891 34221 36925 38654 40743 43368 46182 50036 54043 59324 63100 67446]; % detected infected
D0 = [0 0 0 0 0 1 4 7 11 18 25 34 46 59 77 93 114 136 165 202 244 327 365 445 487 566 688 820 957 1074 1141 1230 1355 1552 1760 1952 2171 2372 2462 2587 2761 2924 3331 3704 4057 4286 4603]; % detected deaths
H0 = [0 1 1 1 2 2 2 2 2 2 2 2 2 2 6 6 6 6 127 127 127 127 127 127 127 127 127 127 173 173 173 173 2979 14026 14026 14026 14026 14026 22130 22991 24325 25318 26573 27655 29160 30152 31142]; % detected recovered

N = 212559417; % population in BR
J = length(I0);
T1 = 4; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % No Imperial
    alpha1 = 5; alpha2 = alpha1;
    pQ = 2; pN = 12; pR = 6;
else
    alpha1 = 12.8; alpha2 = alpha1;
    delay = 20; delay_pc = 14;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(2*delta,2*delta);
input('Press ENTER to proceed, please')
close all;

%% Russia
Country = 'Russia';

% Statistics from 12/03 (Wikipedia + hopkins)
I0 = [45 59 63 93 114 147 199 253 306 367 438 495 658 840 1036 1264 1534 1836 2337 2777 3548 4149 4731 4874 5410 6343 7497 8672 10131 11917 13584 15770 18328 21102 24490 27938 32007 36793 42853 47121 52763 57999 62773 68622 74588 80949 87147 93558]; % detected infected
D0 = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 8 9 9 9 17 24 30 34 43 45 47 58 63 76 94 106 130 148 170 198 232 273 313 361 405 456 513 555 615 681 747 794 876]; % detected deaths
H0 = [3 3 3 4 5 5 5 12 16 16 17 22 29 38 45 49 64 66 121 190 235 281 333 355 395 406 494 580 698 795 1045 1291 1470 1694 1986 2304 2590 3057 3291 3446 3873 4420 4891 5568 6250 6767 7346 8456]; % detected recovered

N = 146745098; % population in RU
J = length(I0);
T1 = 15; % the date of quarantine
T0 = datetime(2020,3,12); % the initial date for data

if (Imperial == 1) % No Imperial
    alpha1 = 5; alpha2 = alpha1;
    pQ = 2; pN = 12; pR = 6;
else
    alpha1 = 2.7; alpha2 = alpha1;
    delay = 20; delay_pc = 20;
end

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(delta,delta);
input('Press ENTER to proceed, please')
close all;

%% China
Country = 'China';

% Statistics from 16/01 (Wikipedia + hopkins)
I0 = [45 62 121 198 291 440 571 830 1287 1975 2744 4515 5974 7711 9692 11791 14380 17205 20438 24324 28018 31161 34546 37198 40171 42638 44653 39571 40364 40381 40007 70548 72436 74185 75002 75891 76288 76936 77150 77658 78064 78497 78824 79251 79824 80026 80151 80270 80409 80552 80651 80695 80735 80754 80778 80793 80813 80824 80844 80860 80881 80894 80928 80967 81008 81054 81093 81171 81218 81285 81340 81394 81439 81470 81518 81554 81589 81620 81639 81669 81708 81740 81802 81865 82924 83014 83134 83135 83302 83351 83402 83753 83784 83804 83817 83849 83864 83876 83884 83899 83909 83912 83938]; % detected infected
D0 = [2 2 3 4 6 9 17 25 41 56 80 106 132 170 213 259 304 361 425 490 563 637 722 811 908 1016 1113 1259 1380 1523 1665 1770 1868 2004 2118 2236 2345 2443 2592 2663 2715 2744 2788 2835 2870 2912 2943 2981 3012 3042 3070 3097 3119 3136 3158 3169 3176 3189 3199 3213 3226 3237 3245 3248 3255 3261 3270 3277 3281 3287 3292 3295 3300 3304 3305 3312 3318 3322 3326 3329 3331 3331 3333 3335 3340 3343 3343 3343 3345 3346 3346 4636 4636 4636 4636 4636 4636 4636 4636 4636 4636 4637 4637]; % detected deaths
H0 = [15 19 24 25 25 25 25 34 38 49 51 60 103 124 171 243 328 475 632 892 1153 1540 2050 2649 3281 3996 4740 5642 6723 8096 9419 10844 12545 14376 16155 18264 20659 22401 24734 27323 29745 32945 36117 39002 41625 44462 47204 49856 52045 53726 55404 57065 58600 59897 61475 62793 64111 65541 66911 67749 68679 69601 70420 71150 71740 72244 72703 73159 73650 74051 74588 74971 75448 75770 76052 76238 76408 76571 76751 76964 77078 77167 77279 77370 77758 77874 77953 77956 78148 78265 78370 77480 77578 77647 77706 77775 77827 77911 77997 78109 78185 78277 78377]; % detected recovered

N = 1438070898; % population in CN https://www.worldometers.info/world-population/china-population/
J = length(I0);
T1 = 13; % the date of quarantine 29/01
T0 = datetime(2020,1,16); % the initial date for data
J_id = 75; %70

alpha1 = 3; alpha2 = alpha1;
delay = 2; delay_pc = 15;
% Confinement parameters for China
pQ = 2; % number of contacts in quarantine
pN = 12; % number of contacts in normal state 
pR = 1; % number of contacts in relaxed quarantine
% Scenario for China
wC = 9; % weeks in quarantine
wR = 6; % weeks in relaxed quarantine
NC = 1; % number cycles to simulate

disp(' ');
disp([' ***' Country '*** ']);
disp(' ');

script(2*delta,2*delta);
disp('Stop execution')