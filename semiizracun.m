## text/x-matlab;
1; % this is script, not a function file
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

function mehkomagnetik = racunajKMM
  global dataL1 dataL2 dataMM
  global KMM KMM1 KMM2

  KMM1 = calcKMM(dataMM,dataL1);
  KMM2 = calcKMM(dataMM,dataL2);
  KMM = KMM1 + KMM2;
  mehkomagnetik = "OK";
endfunction

function polje = racunajPolje( xMin, xMax, zMin, zMax, d )
  global dataL1 dataL2 dataMM
  global KMM KMM1 KMM2
  nTockX  = ceil((xMax-xMin)/d);
  nTockZ  = ceil((zMax-zMin)/d);
  [X, Z]  = meshgrid(
    linspace(xMin,xMax,nTockX),
    linspace(zMin,zMax,nTockZ)
  );
  printf (" - izracun polja v mreži (A, B) [%dx%d tock]\n", nTockZ, nTockX);

  mrezaXZ           = zeros(nTockZ,nTockX);

  A0L1x = mrezaXZ;  B0L1x = mrezaXZ;  B0L1z = mrezaXZ;
  A1L1x = mrezaXZ;  B1L1x = mrezaXZ;  B1L1z = mrezaXZ;
  AL1x  = mrezaXZ;  BL1x  = mrezaXZ;  BL1z  = mrezaXZ;
  A0L2x = mrezaXZ;  B0L2x = mrezaXZ;  B0L2z = mrezaXZ;
  A1L2x = mrezaXZ;  B1L2x = mrezaXZ;  B1L2z = mrezaXZ;
  AL2x  = mrezaXZ;  BL2x  = mrezaXZ;  BL2z  = mrezaXZ;
  A0x   = mrezaXZ;  B0x   = mrezaXZ;  B0z   = mrezaXZ;
  A1x   = mrezaXZ;  B1x   = mrezaXZ;  B1z   = mrezaXZ;
  Ax    = mrezaXZ;  Bx    = mrezaXZ;  Bz    = mrezaXZ;

  for xx = 1:nTockX
  for zz = 1:nTockZ
    T             = [X(zz,xx) 0 Z(zz,xx)];  % Točka v kateri računamo
    % Vektorski magnetni potencial (A)
    A0L1          = A0L(T, dataL1);         % A tuljave 1, brez MM
    A1L1          = A1MM(T, KMM1);          % A vpliva tuljave 1 na MM
    AL1           = A0L1 + A1L1;            % A tuljave 1 z vplivom MM

    A0L2          = A0L(T, dataL2);         % A tuljave 2, brez MM
    A1L2          = A1MM(T, KMM2);          % A vpliva tuljave 2 na MM
    AL2           = A0L2 + A1L2;            % A tuljave 2 z vplivom MM

    A0            = A0L1 + A0L2;            % A obeh tuljav, brez MM
    A1            = A1L1 + A1L2;            % A vplivov obeh tuljav na MM
    A             = A0 + A1;                % A obeh tuljav z vplivom MM

    A0L1x(zz,xx)  = A0L1;   A1L1x(zz,xx)  = A1L1;   AL1x(zz,xx)   = AL1;
    A0L2x(zz,xx)  = A0L2;   A1L2x(zz,xx)  = A1L2;   AL2x(zz,xx)   = AL2;
    A0x(zz,xx)    = A0;     A1x(zz,xx)    = A1;     Ax(zz,xx)     = A;

    % Gostota magnetnega potenciala (B)
    B0L1          = B0L(T, dataL1);         % B tuljave 1, brez MM
    B1L1          = B1MM(T, KMM1);          % B vpliva tuljave 1 na MM
    BL1           = B0L1 + B1L1;

    B0L2          = B0L(T, dataL2);         % B tuljave 2, brez MM
    B1L2          = B1MM(T, KMM2);          % B vpliva tuljave 2 na MM
    BL2           = B0L2 + B1L2;

    B0            = B0L1 + B0L2;            % B obeh tuljav, brez MM
    B1            = B1L1 + B1L2;            % B vplivov obeh tuljav na MM
    B             = B0 + B1;                % B obeh tuljav z vplivom MM
 
    B0L1x(zz,xx)  = B0L1(1);    B0L1z(zz,xx)  = B0L1(3);
    B1L1x(zz,xx)  = B1L1(1);    B1L1z(zz,xx)  = B1L1(3);
    BL1x(zz,xx)   = BL1(1);     BL1z(zz,xx)   = BL1(3);

    B0L2x(zz,xx)  = B0L2(1);    B0L2z(zz,xx)  = B0L2(3);
    B1L2x(zz,xx)  = B1L2(1);    B1L2z(zz,xx)  = B1L2(3);
    BL2x(zz,xx)   = BL2(1);     BL2z(zz,xx)   = BL2(3);

    B0x(zz,xx)    = B0(1);      B0z(zz,xx)    = B0(3);
    B1x(zz,xx)    = B1(1);      B1z(zz,xx)    = B1(3);
    Bx(zz,xx)     = B(1);       Bz(zz,xx)     = B(3);
  endfor
  endfor
  save -z dataPolje X Z A0L1x A1L1x AL1x A0L2x A1L2x AL2x A0x A1x Ax B0L1x B0L1z B1L1x B1L1z BL1x BL1z B0L2x B0L2z B1L2x B1L2z BL2x BL2z B0x B0z B1x B1z Bx Bz
  polje = "OK";
endfunction