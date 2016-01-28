## Copyright (C) 2015 Mitja Jež
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## Author: Mitja Jež <mitja@jež.si>
## Created: 2015-12-14
## Last edit: 2016-01-15
testmode = 0;

if(testmode)
  clear
  testmode = 1;
endif

global mi0 ex ey ez
mi0 = 4*pi*1e-7;            % Permeabilnost praznega prostora
ex  = [1; 0; 0];    ey  = [0; 1; 0];    ez  = [0; 0; 1];
efi = ey;

functions;                  % klic knjiznjice funkcij
semiizracun;

rhoCu   = 8920;   %kg/m³
% ----- računska gostota -------------------------------------------
if(testmode)
  dmin    = 9e-4;   %TESTING
  dmreza  = dmin;   % gostota mreže
  dmag    = dmin;   % debelina segmenta na MM
else
  dmin    = 29e-5;
  dmreza  = 23e-5;  % gostota mreže
  dmag    = 27e-5;  % debelina segmenta na MM
endif
  Sz      = 0.8;    % mm² presek žice Sz = {0.8, 1.5, 2.5, 4}
  dz = 2 * sqrt(Sz*1e-6 / pi); % debelina zice

% ----- računsko območje -------------------------------------------
xMin    = -5e-4;     xMax = 17e-3;
zMin    = -4e-3;     zMax = 20e-3;

% ----- goemetrija in postavitev obeh tuljav --------------
global dataL1 dataL2 dataMM
Nr1     = 1;                % St. plasti ovojev
Nz1     = 13;               % St. ovojev
I1      = -1.3e-3;          % Tok zice v tuljavi
Rin1    = 2e-3;             % Radij tuljavnika
Zoff1   = 11e-3;             % Postavitev tuljave
l1      = Nz1 * dz;         % Dolzina navitja
Rout1   = Rin1 + Nr1*dz;
dataL1  = [I1, l1, Rout1, Rin1, Zoff1, dz];
%printf ("L1: %dx%d ovojev\n", Nz1, Nr1);


Nr2     = 13;               % St. plasti ovojev
Nz2     = 1;                % St. ovojev
I2      = 1.7e-3;           % Tok zice v tuljavi
Rin2    = 3e-3;             % Radij tuljavnika
l2      = Nz2 * dz;         % Dolzina navitja
Zoff2   = l2/2;             % Postavitev tuljave
Rout2   = Rin2 + Nr2*dz;    
dataL2  = [I2, l2, Rout2, Rin2, Zoff2, dz];
%printf ("L2: %dx%d ovojev\n", Nz2, Nr2);

mi      = 1e5;              % permeabilnost MM
lMM     = 17e-3;            % dolzina MM
RMM     = 2e-3;             % radij MM
ZoffMM  = lMM/2;            % položaj MM
NsMM    = ceil(lMM/dmag);   % stevilo segmentov
dlMM    = lMM/NsMM;         % dolzina segmenta
NsOMM   = ceil(RMM/dmag);   % stevilo segmentov
dlOMM   = RMM/NsOMM;        % dolzina segmenta
dataMM  = [mi, lMM, RMM, ZoffMM, NsMM, dlMM, NsOMM, dlOMM];

NsMM    = (NsMM+2*NsOMM);
printf ("~MM: %dx segmentov\n", NsMM);

if(1)
%##########################################################
disp "##### izracunavanje ###############################";
%##########################################################

polje = 0;
sile = 1;
induktivnosti = 0;

if(exist ("mehkomagnetik", "var") == 0 )
  global KMM KMM1 KMM2
  printf (" - mehkomagnetik (KMM)\n");
  [KMM_start T1] = timestring();
  mehkomagnetik = racunajKMM;
  [KMM_end T2 hhmmss] = timestring(T1);
  printf ("   (izračunal v času: %1.0f h, %2.0f min in %2.0f s)\n", hhmmss(1), hhmmss(2), hhmmss(3));
else
  printf(" - tokovno polje na MM je ŽE %s\n", mehkomagnetik);
endif %mreza


if( induktivnosti)
  if( exist ("Wm", "var") == 0 )
    [energija_start T1]  = timestring();
    printf (" - emergija (");
    Wm = getWm(KMM);
    printf ("Wm = %1.2e", Wm);
    printf (")\n");
    [energija_end T2 hhmmss] = timestring(T1);
    printf ("   (izračunal v času: %d h, %d min in %d s)\n", hhmmss(1), hhmmss(2), hhmmss(3));
  endif
endif

if( induktivnosti )
  if( exist ("L12", "var") == 0 )
    printf (" - induktivnosti (");
    [L_start T1] = timestring();
    L1 = lastnaL(dataL1, KMM1);
    printf ("L11 = %1.2e H", L1);
    printf (", ");
    L2 = lastnaL(dataL2, KMM2);
    printf ("L22 = %1.2e H", L2);
    printf (", ");
    L12 = (Wm - L1*I1^2 - L2*I2^2)/(I1*I2);
    printf ("L12 = %1.2e H", L12);
    printf (")\n");
    [L_end T2 hhmmss] = timestring(T1);
    printf ("   (izračunal v času: %d h, %d min in %d s)\n", hhmmss(1), hhmmss(2), hhmmss(3));
  endif
endif

if(polje && exist ("izracunanoPolje", "var") == 0 )
  [mreza_start T1]  = timestring();
  izracunanoPolje = racunajPolje( xMin, xMax, zMin, zMax, dmreza )
  [mreza_end T2 hhmmss] = timestring(T1);
  printf ("   (izračunal v času: %1.0f h, %2.0f min in %2.0f s)\n", hhmmss(1), hhmmss(2), hhmmss(3));
elseif (polje)
  printf (" - polje v mreži je ŽE %s\n", izracunanoPolje);
endif %mreza




if( sile && exist("izracunaneSile", "var") == 0 )
  global KMM
  disp " - Wm v različnih zamikih";

  dlF = 5e-4;                             % zamik tuljave L1
  
  [sila_start T1]       = timestring();
  global dataL1 dataL2 dataMM
  lMM     = dataMM(2);                    % končna vrenosto zamika je enaka dolžini MM
  lL      = dataL1(2);
  Z0      = -1.5*lL:dlF:lMM+1.5*lL;       % serija zamikov tuljave
  NF      = numel(Z0);                    % Št. vrednosti
  dataL1tmp = dataL1;
  for n = 1:NF;
    dataL1(5) = Z0(n);                    % Zamik tuljave - Zoff
    ok        = racunajKMM;
    Wm        = getWm(KMM);
    ArrWm(n)  = Wm;
%  printf ("Zamik %d: Wm = %1.2e J \n", n, Wm);
  endfor
  dataL1 = dataL1tmp;
  Fz = diff(ArrWm)/dlF;
  [sila_end T2 hhmmss]  = timestring(T1);
  printf ("   (izračunal v času: %d h, %d min in %d s)\n", hhmmss(1), hhmmss(2), hhmmss(3));
  izracunaneSile = "OK";

% SHRANI IZRAČUNE
  save -z dataSile dlF Z0 ArrWm Fz
endif %izracunSile





if(1)
  semigraf %izris
endif %izris

endif %izracunavanje
