% Site data: Mf prevalence, blood sample volume, MDA regimen, MDA frequency,...
% MDA coverage, number of years of treatment, vector control, switch year

function [ DokanTofaMf,PiapungMf,SeriMf,GbuwhenMf,MaigangaMf,KirareMf,GizaMf,UsinoBundiMf,...
    DokanTofaCfa,PiapungCfa,SeriCfa,GbuwhenCfa,MaigangaCfa,KirareCfa,GizaCfa,UsinoBundiCfa,...
    DokanTofaL3,PiapungL3,SeriL3,GbuwhenL3,MaigangaL3,KirareL3,GizaL3,UsinoBundiL3,...
    DokanTofaMfAge,PiapungMfAge,SeriMfAge,GbuwhenMfAge,MaigangaMfAge,KirareMfAge,GizaMfAge,UsinoBundiMfAge,...
    DokanTofaCfaAge,PiapungCfaAge,SeriCfaAge,GbuwhenCfaAge,MaigangaCfaAge,KirareCfaAge,GizaCfaAge,UsinoBundiCfaAge,...
    DokanTofaYear,PiapungYear,SeriYear,GbuwhenYear,MaigangaYear,KirareYear,GizaYear,UsinoBundiYear,...
    DokanTofaReg,PiapungReg,SeriReg,GbuwhenReg,MaigangaReg,KirareReg,GizaReg,UsinoBundiReg,...
    DokanTofaFreq,PiapungFreq,SeriFreq,GbuwhenFreq,MaigangaFreq,KirareFreq,GizaFreq,UsinoBundiFreq,...
    DokanTofaMDACov,PiapungMDACov,SeriMDACov,GbuwhenMDACov,MaigangaMDACov,KirareMDACov,GizaMDACov,UsinoBundiMDACov,...
    DokanTofaNumYears,PiapungNumYears,SeriNumYears,GbuwhenNumYears,MaigangaNumYears,KirareNumYears,GizaNumYears,UsinoBundiNumYears,...
    DokanTofaVC,PiapungVC,SeriVC,GbuwhenVC,MaigangaVC,KirareVC,GizaVC,UsinoBundiVC,...
    DokanTofaSwitchYear,PiapungSwitchYear,SeriSwitchYear,GbuwhenSwitchYear,MaigangaSwitchYear,KirareSwitchYear,GizaSwitchYear,UsinoBundiSwitchYear,...
    DokanTofaITNCov,PiapungITNCov,SeriITNCov,GbuwhenITNCov,MaigangaITNCov,KirareITNCov,GizaITNCov,UsinoBundiITNCov,...
    DokanTofaITNYear,PiapungITNYear,SeriITNYear,GbuwhenITNYear,MaigangaITNYear,KirareITNYear,GizaITNYear,UsinoBundiITNYear,...
    DokanTofaIRSCov,PiapungIRSCov,SeriIRSCov,GbuwhenIRSCov,MaigangaIRSCov,KirareIRSCov,GizaIRSCov,UsinoBundiIRSCov]...
    = PostIntv_data_SSA
%% Post-intervention Mf and Cfa Data
% enter overall community prevalence in a single line
% 1st column = month of survey; 2nd: Total number of samples; 3rd: Mf +ves;

%%% Mf
DokanTofaMf = [(2005-2002)*12+1 236 7
    (2006-2002)*12+1 132 0
    (2007-2002)*12+1 151 2
    (2008-2002)*12+1 158 0
    (2009-2002)*12+1 223 1
    (2010-2002)*12+1 119 0
    (2011-2002)*12+1 206 1];
PiapungMf = [(2005-2002)*12+1 256 11
    (2006-2002)*12+1 160 6
    (2007-2002)*12+1 187 18
    (2009-2002)*12+1 291 6
    (2011-2002)*12+1 90 1];
SeriMf = [%(2003-1995)*12+1 527 56 % age profile
    (2004-1995)*12+1 97 11
    (2005-1995)*12+1 321 5
    (2006-1995)*12+1 157 2
    (2007-1995)*12+1 133 1
    (2008-1995)*12+1 110 3
    (2009-1995)*12+1 258 0
    (2010-1995)*12+1 268 3
    (2011-1995)*12+1 211 0];
GbuwhenMf = [%(2000-1992)*12+1 30 1
    % (2003-1992)*12+1 508 19 % age profile
    (2004-1992)*12+1 446 8
    (2005-1992)*12+1 286 1
    (2006-1992)*12+1 183 1
    (2007-1992)*12+1 197 0
    (2008-1992)*12+1 127 0
    (2009-1992)*12+1 175 0
    (2010-1992)*12+1 127 0
    (2011-1992)*12+1 127 0];
MaigangaMf = [%(2000-1992)*12+1 109 5
    % (2003-1992)*12+1 478 23 % age profile
    (2004-1992)*12+1 286 15
    (2005-1992)*12+1 169 5
    (2006-1992)*12+1 126 7
    (2007-1992)*12+1 158 1
    (2008-1992)*12+1 109 2
    (2009-1992)*12+1 152 1
    (2010-1992)*12+1 109 0
    (2011-1992)*12+1 109 0];
KirareMf = [%(2004-2004)*12+1 471 round(0.261*471);
               (2005-2004)*12+1 461 96;
               (2006-2004)*12+1 438 69;
               (2007-2004)*12+1 351 35;
               (2008-2004)*12+1 302 39;
               (2009-2004)*12+1 259 13;
               (2010-2004)*12+1 400 17;
               (2011-2004)*12+1 393 11;
               (2013-2004)*12+1 60 round(0.055*60)];
PenengMf = [%(1994-1994)*12+1 63 round(0.667*63);
               (1995-1994)*12+1 65 round(0.615*65);
               (1996-1994)*12+1 88 round(0.205*88);
               (1997-1994)*12+1 89 round(0.135*89);
               (1998-1994)*12+1 92 round(0.054*92);
               (1999-1994)*12+1 109 round(0.037*109)];
GizaMf = [%0*12+1 1067 123;
           1*12+1 1012 46;
           2*12+1 1026 28;
           3*12+1 1010 13;
           4*12+1 1116 4;
           5*12+1 1064 13];
UsinoBundiMf = [%0*12+1 757 106;
               1*12+1 696 58;
               2*12+1 714 24;
               3*12+1 529 7];

%%% Cfa
DokanTofaCfa = [(2007-2002)*12+1 277 40
    (2008-2002)*12+1 158 14
    (2009-2002)*12+1 223 7
    (2010-2002)*12+1 119 3
    (2011-2002)*12+1 206 7];
PiapungCfa = [(2005-2002)*12+1 192 46
    (2007-2002)*12+1 312 61
    (2008-2002)*12+1 62 9
    (2009-2002)*12+1 291 27
    (2011-2002)*12+1 90 2];
SeriCfa = [% (2003-1995)*12+1 528 178 % age profile
    (2004-1995)*12+1 474 106
    (2005-1995)*12+1 150 20
    (2007-1995)*12+1 385 87
    (2008-1995)*12+1 110 9
    (2009-1995)*12+1 258 27
    (2010-1995)*12+1 268 34
    (2011-1995)*12+1 211 7
    (2016-1995)*12+1 74 0]; 
GbuwhenCfa = [%(2000-1992)*12+1 30 14
    % (2003-1992)*12+1 510 69 % age profile
    (2004-1992)*12+1 446 23
    (2005-1992)*12+1 178 8
    (2007-1992)*12+1 127 8
    (2008-1992)*12+1 175 0
    (2009-1992)*12+1 127 1
    (2010-1992)*12+1 127 1
    (2011-1992)*12+1 127 1];
MaigangaCfa = [%(2000-1992)*12+1 50 27
    % (2003-1992)*12+1 487 87 % age profile
    (2004-1992)*12+1 50 8
    (2005-1992)*12+1 165 34
    (2007-1992)*12+1 109 21
    (2008-1992)*12+1 152 3
    (2009-1992)*12+1 50 4
    (2010-1992)*12+1 50 6
    (2011-1992)*12+1 50 2];
KirareCfa = [%(2004-2004)*12+1 90 round(0.261*471);
               (2005-2004)*12+1 66 35;
               (2006-2004)*12+1 72 37;
               (2007-2004)*12+1 61 32;
               (2008-2004)*12+1 49 22;
               (2010-2004)*12+1 400 101;
               (2011-2004)*12+1 393 77;
               (2013-2004)*12+1 422 69];
PenengCfa = NaN;
GizaCfa = [%0*12+1 1067 203;
           1*12+1 1012 162;
           2*12+1 1026 108;
           3*12+1 1010 53;
           4*12+1 1116 31;
           5*12+1 1064 51];
UsinoBundiCfa = [%0*12+1 558 265;
               1*12+1 692 243;
               2*12+1 695 175;
               3*12+1 543 93];

DokanTofaL3 = [(2003-2002)*12+1 1241	12
(2004-2002)*12+1 741	6
(2005-2002)*12+1 897	51
(2006-2002)*12+1 548	26
(2007-2002)*12+1 402	5
(2008-2002)*12+1 346	4
(2009-2002)*12+1 684	1
(2010-2002)*12+1 507	1
(2011-2002)*12+1 176	0
(2012-2002)*12+1 119	0
(2013-2002)*12+1 102	0];

PiapungL3 = [(2003-2002)*12+1 2	0
(2004-2002)*12+1 688	49
(2005-2002)*12+1 718	15
(2006-2002)*12+1 768	77
(2007-2002)*12+1 547	3
(2008-2002)*12+1 306	2
(2009-2002)*12+1 479	1
(2010-2002)*12+1 183	0
(2011-2002)*12+1 130	0
(2012-2002)*12+1 194	0
(2013-2002)*12+1 104	0];

SeriL3 = [(2000-1995)*12+1 327	13
(2001-1995)*12+1 373	42
(2002-1995)*12+1 184	4
(2003-1995)*12+1 975	24
(2004-1995)*12+1 479	14
(2005-1995)*12+1 438	7
(2006-1995)*12+1 379	2
(2007-1995)*12+1 495	396
(2008-1995)*12+1 389	6
(2009-1995)*12+1 579	6
(2010-1995)*12+1 276	3
(2011-1995)*12+1 225	0
(2012-1995)*12+1 194	0
(2013-1995)*12+1 163	0];

GbuwhenL3 = [(2000-1992)*12+1 255	2
(2001-1992)*12+1 460	2
(2002-1992)*12+1 391	0
(2003-1992)*12+1 1840	6
(2004-1992)*12+1 1285	7
(2005-1992)*12+1 1119	1
(2006-1992)*12+1 858	2
(2007-1992)*12+1 732	2
(2008-1992)*12+1 505	1
(2009-1992)*12+1 569	0
(2010-1992)*12+1 566	0
(2011-1992)*12+1 407	0];

MaigangaL3 = [(2000-1992)*12+1 73	10
(2001-1992)*12+1 378	12
(2002-1992)*12+1 397	3
(2003-1992)*12+1 2202	10
(2004-1992)*12+1 842	8
(2005-1992)*12+1 923	8
(2006-1992)*12+1 976	5
(2007-1992)*12+1 568	1
(2008-1992)*12+1 583	0
(2009-1992)*12+1 701	2
(2010-1992)*12+1 447	0
(2011-1992)*12+1 120	0];

KirareL3 = [(2004-2004)*12+1 5396	77
(2005-2004)*12+1 7840	35
(2006-2004)*12+1 7071	12
(2007-2004)*12+1 9088	10	
(2009-2004)*12+1 1830	1
(2010-2004)*12+1 3976	4
(2011-2004)*12+1 7713	1
(2012-2004)*12+1 125	0
(2013-2004)*12+1 1532	0];

GizaL3 = [];

UsinoBundiL3 = [];

%% Post-intervention Mf age profile
DokanTofaMfAge = NaN;  
PiapungMfAge = NaN;  
SeriMfAge = [4.5 91 1 9
    14.5 194 14 19
    24.5 102 15 29
    34.5 26 4 39
    44.5 78 15 49
    54.5 28 4 59
    64.5 9 3 69];
GbuwhenMfAge = [6 198 2 10  
    14.5 99 1 20
    24.5 86 6 30
    34.5 57 4 40	
    44.5 27 2 50	
    54.5 16 4 60	
    64.5 8 0 70];   
MaigangaMfAge = [6 174 1 10  
    14.5 111 5 20
    24.5 69 7 30
    34.5 49 2 40	
    44.5 26 4 50	
    54.5 19 2 60	
    64.5 21 2 70];  
KirareMfAge = NaN;
PenengMfAge = NaN;
GizaMfAge = NaN;
UsinoBundiMfAge = NaN;

DokanTofaCfaAge = NaN;  
PiapungCfaAge= NaN;  
SeriCfaAge = [4.5 91 4 9
    14.5 194 58 19
    24.5 102 44 29
    34.5 26 13 39
    44.5 78 41 49
    54.5 28 14 59
    64.5 9 4 69];
GbuwhenCfaAge = [6 198 9 10  
    14.5 101 12 20
    24.5 86 20 30
    34.5 57 13 40	
    44.5 27 5 50	
    54.5 16 6 60	
    64.5 8 4 70]; 
MaigangaCfaAge = [6 174 1 10  
    14.5 111 5 20
    24.5 69 7 30
    34.5 49 2 40	
    44.5 26 4 50	
    54.5 19 2 60	
    64.5 21 2 70];
KirareCfaAge = NaN;
PenengCfaAge = NaN;
GizaCfaAge = NaN;
UsinoBundiCfaAge = NaN;

% year of first age profile
DokanTofaYear = (2007-2002)+1; % year to create age profile for hindcasting
PiapungYear = (2005-2002)+1; % year to create age profile for hindcasting
SeriYear = (2003-1995)+1;
GbuwhenYear = (2003-1992)+1;
MaigangaYear = (2003-1992)+1;
KirareYear = 4; % year to create age profile for hindcasting
PenengYear = 1;
GizaYear = 3; % year to create age profile for hindcasting
UsinoBundiYear = 2; % year to create age profile for hindcasting

%% MDA Regimen
% if more than one treatment regimen, given frequency for each regimen
% (i.e. [5,2] for ALB followed by DEC+ALB)

% 1: IVM
% 2: DEC+ALB
% 3: IVM+DEC+ALB
% 4: IVM+ALB
% 5: ALB
% 6: DEC salt

DokanTofaReg = [4];
PiapungReg = [4];
SeriReg = [1,4];
GbuwhenReg = [1,4];
MaigangaReg = [1,4];
KirareReg = [4];
PenengReg = [7];
GizaReg = [2];
UsinoBundiReg = [2];

%% MDA Frequency
% in months
% if more than one treatment regimen, given frequency for each regimen
% (i.e. [12,6] for annual followed by biannial)

DokanTofaFreq = [12];
PiapungFreq = [12];
SeriFreq = [12,12];
GbuwhenFreq = [12,12];
MaigangaFreq = [12,12];
KirareFreq = [12];
PenengFreq = [12];
GizaFreq = [12];
UsinoBundiFreq = [12];

%% Annual MDA Coverage
% enter as proportion (i.e. 0.8)

% MDA zeros included in order to model up to latest data points
DokanTofaMDACov = [0 74.9 76.7 67.4 77.6 77.1 78.3 78.2 0 0 0 0]/100; % 2002-2013 
PiapungMDACov = [0 70.2 72.0 78.0 78.5 80.1 79.2 78.9 0 0 0 0]/100; % 2002-2013
SeriMDACov = [77.1 77.7	79.3 79.1 68.4 75.2 74.1 76.7 78.2 79.3 78.1 ...
    75.9 87.4 75.0 85.1 91.1 81.1 78.7 79.7 79.7 79.8 79.7 0]/100; % 1995-2017
GbuwhenMDACov = [52.4 47.8 41.8 67.5 89.1 85.8 51.9 62.3 68.0 62.7 49.6...
    98.7 93.9 96.6 92.5 91.8 93.9 79.3 88.3 88.3 88.3 88.3 88.3 88.3 88.3 0]/100; % 1992-2017
MaigangaMDACov = [70.3 85.8 85.6 89.5 94.6 77.4 84.6 73.6 74.1 85.3 81.6...
    79.8 82.8 90.5 91.5 92.8 72.9 89.7 85.1 85.1 85.1 85.1 85.1 85.1 85.1 0]/100; % 1992-2017
KirareMDACov = [0.64 0.76 0.696 0 0.773 0.789 0.6 0.408 0.368 0];
% PenengMDACov = [0.50 0.78 0.75 0.68 0.72];
GizaMDACov = [0.867 0.955 0.901 0.888 0.903 0];
UsinoBundiMDACov = [0.684 0.764 0.739 0 ];


%% Total number of years of treatment and/or data
% (last year - first year) + 1
DokanTofaNumYears = (2013-2002)+1;
PiapungNumYears = (2013-2002)+1;
SeriNumYears = (2017-1995)+1;
GbuwhenNumYears = (2017-1992)+1;
MaigangaNumYears = (2017-1992)+1;
KirareNumYears = (2013-2004)+1;
% PenengNumYears = (1999-1994)+1;
GizaNumYears = 6;
UsinoBundiNumYears = 4;



%% Vector control
% 0: no VC, 1: VC

DokanTofaVC = 1;
PiapungVC = 1;
SeriVC = 1;
GbuwhenVC = 1;
MaigangaVC = 1;
KirareVC = 1;
% PenengVC = 0;
GizaVC = 0;
UsinoBundiVC = 1;

%% Year to switch treatment regimens
% enter 0 if a single treatment regimen is followed

DokanTofaSwitchYear = 0;
PiapungSwitchYear = 0;
SeriSwitchYear = 6;
GbuwhenSwitchYear = 8;
MaigangaSwitchYear = 9;
KirareSwitchYear = 0;
% PenengSwitchYear = 0;
GizaSwitchYear = 0;
UsinoBundiSwitchYear = 0;

%% Vector control Coverage
% enter as proportion (i.e. 0.8)
DokanTofaITNYear = (2004-2002)+1;
PiapungITNYear = (2004-2002)+1;
SeriITNYear = (2004-1995)+1;
GbuwhenITNYear = (2004-1992)+1;
MaigangaITNYear = (2004-1992)+1;
KirareITNYear = 1;
GizaITNYear = 1;
UsinoBundiITNYear = 1;

DokanTofaITNCov = 0.5; 
PiapungITNCov = 0.5; 
SeriITNCov = 0.804;
GbuwhenITNCov = 0.5; 
MaigangaITNCov = 0.5; 
KirareITNCov = 0.25;
% PenengITNCov = 0;
GizaITNCov = 0;
UsinoBundiITNCov = 0.4; 

DokanTofaIRSCov = 0;
PiapungIRSCov = 0;
SeriIRSCov = 0;
GbuwhenIRSCov = 0;
MaigangaIRSCov = 0;
KirareIRSCov = 0;
% PenengIRSCov = 0;
GizaIRSCov = 0;
UsinoBundiIRSCov = 0;

end
