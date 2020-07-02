function [ DokanTofaMf,PiapungMf,SeriMf,GbuwhenMf,MaigangaMf,KirareMf,PenengMf,GizaMf,UsinoBundiMf,...
    DokanTofaCfa,PiapungCfa,SeriCfa,GbuwhenCfa,MaigangaCfa,KirareCfa,PenengCfa,GizaCfa,UsinoBundiCfa,...
    DokanTofaVol,PiapungVol,SeriVol,GbuwhenVol,MaigangaVol,KirareVol,PenengVol,GizaVol,UsinoBundiVol,...
    DokanTofaABR,PiapungABR,SeriABR,GbuwhenABR,MaigangaABR,KirareABR,PenengABR,GizaABR,UsinoBundiABR,...
    DokanTofabCulex,PiapungbCulex,SeribCulex,GbuwhenbCulex,MaigangabCulex,KirarebCulex,PenengbCulex,GizabCulex,UsinoBundibCulex]...
    = baseline_data

%% Baseline Mf Data
% Best data to provide would be age-stratified prevalence
% Variable names = {SiteName}Mf

%% pre-intervention mf and cfa data
DokanTofaMf = [4.5 49 1 9
    14.5 137 4 19
    24.5 101 8 29
    34.5 61	4 39
    44.5 33	2 49
    54.5 24	1 59
    64.5 10	1 69
    74.5 4	0 79];
PiapungMf = [4.5 62	1 9
    14.5 100 3 19
    24.5 74	5 29
    34.5 66	8 39
    44.5 43	6 49
    54.5 36	12 59
    64.5 14	3 69
    74.5 8	2 79];
SeriMf = NaN;
GbuwhenMf = NaN;
MaigangaMf = NaN;
KirareMf = [4	281	23	9
    14.5	206	58	19
    24.5	111	42	29
    34.5	115	37	39
    44.5	74	21	49
    59.5	132	44	69];
PenengMf = [5	7	2   10
    15	14	9   20
    25	18	11  30
    35	10	8   40
    45	8	8   50
    55	6	4   60];
GizaMf = [35 1067 123 70];
UsinoBundiMf = [8 125 3 10
    13 90 7 15
    18 60 13  20
    25.5 131 31 30
    35.5 80 24 40
    45.5 49 24 50
    60.5 37 13 70];


DokanTofaCfa = [4.5 49 2 9
    14.5 137 14 19
    24.5 101 27 29
    34.5 61	32 39
    44.5 33	6 49
    54.5 24	11 59
    64.5 10	4 69
    74.5 4 1 79];
PiapungCfa = [4.5 62 3 9
    14.5 100 13 19
    24.5 74	30 29
    34.5 66	23 39
    44.5 43	18 49
    54.5 36	25 59
    64.5 14	7 69
    74.5 8 3 79];
SeriCfa = NaN;
GbuwhenCfa = NaN;
MaigangaCfa = NaN;
KirareCfa = [34.5 90 48 69];
PenengCfa = NaN;
GizaCfa = [35 1067 203 70];
UsinoBundiCfa = [8 120 22 10
    13 87 25 15
    18 60 27  20
    25.5 128 78 30
    35.5 79 52 40
    45.5 47 34 50
    60.5 37 27 70];

%% Blood sample volume used to test presence of mf (in uL)
% Standard volumes are 1000 uL, 100 uL, or 20 uL (60 uL uses 100 uL
% correction)
% Model operates based on 1 mL samples, so smaller samples will be
% corrected to reflect this

DokanTofaVol = 60;
PiapungVol = 60;
SeriVol = 60;
GbuwhenVol = 60;
MaigangaVol = 60;
KirareVol = 100;
PenengVol = 1000;
GizaVol = 1000;
UsinoBundiVol = 1000;

%% Baseline ABR data
% If not available for a particular site, enter value as NaN

DokanTofaABR = NaN;  %% simplified into one column
PiapungABR = NaN; 
SeriABR = NaN; 
GbuwhenABR = NaN;
MaigangaABR = NaN;
KirareABR = 2091;
PenengABR = 8194;
GizaABR = NaN;
UsinoBundiABR = NaN;


%% Mosquito species flag based on dominant genera (0: Anopheles, 1: Culex)

DokanTofabCulex = 0;  %% simplified into one colume,  name change to specy.
PiapungbCulex = 0; 
SeribCulex = 0;
GbuwhenbCulex = 0;
MaigangabCulex = 0;
KirarebCulex = 0;
PenengbCulex = 0;
GizabCulex = 1;
UsinoBundibCulex = 0;


end
