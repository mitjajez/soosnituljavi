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

functions;
%##########################################################
disp "##### izris #######################################";                                             #
%##########################################################
[izris_start T1] = timestring();

%close;           % zapre zadnje odprto okno
close all;				% zapre vsa odprta okna

global dataL1 dataL2 dataMM

startx = 0; 			%levo
starty = 0; 			%spodej
width  = 1366;    %px
height = 700;     %px

if(polje)
load dataPolje
if(exist ("izracunanoPolje", "var") == 0)
  disp " polje še ni poračunano"
endif
A0  = 1;
A1  = 0;
A   = 0;
row = A0 + A1 + A;

figure('Position',[startx,starty,width,height]);	% odpremo okno
if(A0)
p = 1;
subplot (row,3,p,'align');
  disp "   ~ contour (A0L1)";
	cA = contour(X,Z,A0L1x,12);
  hold on;
  narisiTuljavo(dataL1,1);
  narisiTuljavo(dataL2);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A0L1: Osnovni vektorski magnetni potencial prve tuljave');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');

p = 2;
subplot (row,3,p,'align');
  disp "   ~ contour (A0L2)";
	cA0 = contour(X,Z,A0L2x,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A0L2: Osnovni vektorski magnetni potencial druge tuljave');
	axis ('equal');	xlabel('{\itx} / m');
  
p = 3;
subplot (row,3,p,'align');
  disp "   ~ contour (A0 )";
	cA0 = contour(X,Z,A0x,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A0: Skupni osnovni vektorski magnetni potencial tuljav');
	axis ('equal');	xlabel('{\itx} / m');

  print -dsvg plotAL1.svg
endif

if(A1)
p = 3*A0 + 1;
subplot (row,3,p,'align');
  disp "   ~ contour (A0L2)";
	cA = contour(X,Z,A1L1x,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2,1);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A1L1: Vektorski magnetni potencial vpliva magnetika zaradi prve tuljave');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');

p = 3*A0 + 2;
subplot (row,3,p,'align');
  disp "   ~ contour (A1L2)";
	cA0 = contour(X,Z,A1L2x,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A1L2: Vektorski magnetni potencial vpliva magnetika zaradi druge tuljave');
	axis ('equal');	xlabel('{\itx} / m');
  
p = 3*A0 + 3;
subplot (row,3,p,'align');
  disp "   ~ contour (A1 )";
	cA0 = contour(X,Z,A1x,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A1: Skupni vektorski magnetni potencial vpliva magnetika zaradi obeh tuljav');
	axis ('equal');	xlabel('{\itx} / m');
endif

if(A)
p = 3*(row-1) + 1;
subplot (row,3,p,'align');
  disp "   ~ contour (AL1)";
	cA = contour(X,Z,AL1x,12);
  hold on;
  narisiTuljavo(dataL1,1);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('AL1: Skupni vektorski magnetni potencial prve tuljave');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');

p = 3*(row-1) + 2;
subplot (row,3,p,'align');
  disp "   ~ contour (AL2)";
	cA0 = contour(X,Z,AL2x,12);
  hold on;
  narisiTuljavo(dataL2,1);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('AL2: Skupni vektorski magnetni potencial druge tuljave');
	axis ('equal');	xlabel('{\itx} / m');

p = 3*(row-1) + 3;
subplot (row,3,p,'align');
  disp "   ~ contour (A)";
	cA0 = contour(X,Z,Ax,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  line([0 0],[zMin zMax],'LineWidth',1,'Color',[0.5  0.5  0.5],'LineStyle','-.');
  hold off;
	title ('A: Skupni vektorski magnetni potencial obeh tuljav');
	axis ('equal');	xlabel('{\itx} / m');
  
endif
if(row)
print -dsvg plotA.svg;
endif

endif
if(sile)
load dataSile
global dataL1 dataL2 dataMM
kF      = 1; %1e7;
N       = numel(ArrWm);
Fz(N)   = 0;
ArrWm   = kF.*ArrWm;
lMM     = dataMM(2);
OffMM   = dataMM(4);
lL2     = dataL2(2);
OffL2   = dataL2(5);
Wmax    = max(ArrWm);
Wmean   = mean(ArrWm);
Wmin    = min(ArrWm);
% ----------------------------------------------------------------------
figure('Position',[startx,starty,width,height]);	% odpremo okno
  disp "   ~ Magnetna energija pri več zamikih (ArrWm)";
  plot (Z0, ArrWm, 'r-');
  hold on;
  qB = quiver(Z0,ArrWm,Fz,0);
  line ([OffMM-lMM/2,OffMM+lMM/2], [Wmean,Wmean],'LineWidth',2,'Color','k','LineStyle','-');
  line ([OffL2-lL2/2,OffL2+lL2/2], [Wmin,Wmin],'LineWidth',2,'Color','k','LineStyle','-');
  hold off;
	title ('ArrWm: Magnetna energija pri več zamikih');
  axis ([0,1,0.8*Wmin,1.2*Wmax]);
  axis ('ticx');
  axis ('autox');
  xlabel('{\itz} / m'); ylabel('{{\itW}_{m}}');
  set (gca, 'xaxislocation', 'zero');
  set (gca, 'yaxislocation', 'zero');
  box off

 
print -dsvg plotArrWm.svg
endif
















if(0)
% ----------------------------------------------------------------------
figure('Position',[startx,starty,width,height]);	% odpremo okno
  disp "   ~ quiver  (B0)";
  qB = quiver(X,Z,B0x,B0z);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  hold off;
	title ('B0: Gostota magnetnega pretoka tuljave - osnovnega polja');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');
endif

if(0)
% ----------------------------------------------------------------------
figure('Position',[startx,starty,width,height]);	% odpremo okno
subplot (1,2,1);
  disp "   ~ contour (A1)";
	cA0 = contour(X,Z,A1x,12);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  for ii = 1:NsMM
    K = KMM(ii);
    [T n] = TnMM(ii);     text (T(1), T(3), sprintf("%d",ii));
    TnMMd(ii,1) = T(1);   TnMMd(ii,2) = T(3);
    TnMMd(ii,3) = n(1)*K;   TnMMd(ii,4) = n(3)*K;
  endfor
  qG = quiver( TnMMd(:,1), TnMMd(:,2), TnMMd(:,3), TnMMd(:,4) );
  plot (TnMMd(:,1), TnMMd(:,2), 'r@x');
  hold off;
  colorbar ('southoutside');
	title ('A1: Vektorski magnetni potencial tuljave - vpliv magnetika');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');

subplot (1,2,2);
  disp "   ~ quiver  (B1)";
  qB = quiver(X,Z,B1x,B1z);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  for ii = 1:NsMM
    K = KMM(ii);
    [T n] = TnMM(ii);     text (T(1), T(3), sprintf("%d",ii));
    TnMMd(ii,1) = T(1);   TnMMd(ii,2) = T(3);
  endfor
  plot (TnMMd(:,1), TnMMd(:,2), 'r@x');
  hold off;
	title ('B1: Gostota magnetnega pretoka magnetikov - odzivnega polja');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');
endif

if(0)
% ----------------------------------------------------------------------
figure('Position',[startx,starty,width,height]);	% odpremo okno
  disp "   ~ quiver  (B )";
  qB  = quiver(X,Z,Bx,Bz);
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  hold off;
	title ('A: Vektorski magnetni potencial tuljave - z vplivom magnetika');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');
endif

if(0)
% ----------------------------------------------------------------------
figure('Position',[startx,starty,width,height]);	% odpremo okno
  disp "   ~ geometrija";
  printf("%d", NsMM);
  for ii = 1:NsMM
    K = KMM(ii);
    [T n] = TnMM(ii);     text (T(1), T(3), sprintf("%d",ii));
    TnMMd(ii,1) = T(1);   TnMMd(ii,2) = T(3);
    TnMMd(ii,3) = n(1)*K;   TnMMd(ii,4) = n(3)*K;
  endfor
%  qG = quiver( TnMMd(:,1), TnMMd(:,2), TnMMd(:,3), TnMMd(:,4) );
  q7 = quiver(vec(:,1),vec(:,2),vec(:,3),vec(:,4));
  hold on;
  narisiTuljavo(dataL1);
  narisiTuljavo(dataL2);
  narisiMM(dataMM);
  plot (TnMMd(:,1), TnMMd(:,2), 'r@x');
  hold off;
	title ('GEOMETRIJA');
	axis ('equal');	xlabel('{\itx} / m'); ylabel('{\itz} / m');
endif %izracunMreze
