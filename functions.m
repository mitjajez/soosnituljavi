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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%! @file function.m
%! @bief Funkcije za izračun osnovnega magnetnega polja B0 posamezne tuljave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R m]  = R0_m(r,z,ro,zo)
  % za porazdelitev vektorja gostote magnetnega pretoka
  % (r, z ) ... splosna tocka
  % (ro,zo) ... ovoj, znaka
  Rn2 = (r+ro)^2 + (z-zo)^2;
  m = 4*r*ro/( Rn2 );
  R = sqrt(Rn2);
endfunction


% ----------------------------------------------------------------------
function G = fG(m,R,k)
%/**
% * Greenova funkcija, ki pripada kroznemu ovoju
% * ELMG str. 201
% */
  [K,E] = ellipke(m);  % Eliptic integrala
  G = (2-m)*K - 2*E;
  G = G/(m*R);
  if(k)
    G = G/(2*pi*pi)
  endif
endfunction

% ----------------------------------------------------------------------
function A0 = A0L(T, data)
%/**
% * Vektorski magnentni potencial v tocki (r,z)
% * ELMG str. 201
% */
  global mi0
  [I, l, Rout, Rin, Zoff, d] = vec2var(data);
  [RR ZZ NN] = osiZank(data);
  [r fi z] = kks2vks(T);
  Rz = d/2;                         % radij zice
  A0 = 0;

  if r > 0                          % A0(r=0) = 0
    for n = 1:NN;                   % po vseh zankah
    Ri = RR(n);   Zi = ZZ(n);
      t = sqrt( (Ri-r)^2 + (Zi-z)^2 ); 
      if(t <= Rz)                   %str 199
        A0 += (log(8*Ri/Rz) - 2) + 0.5*(1 - (t/Rz)^2);
        % k = mi0*I/(2*pi)
      else
        [R0 m]  = R0_m(r,z, Ri,Zi);
        A0 += 2*fG(m,R0,0);         % superpozicija ovojev, fG brez konstant
        % k = 1/(2*pi^2)
        % Afi = 2*pi*mi0*I * fG(m,R0,1)% racunamo enacbo na str.: 201 -----
        % Afi = 2*pi*mi0*I * 1/2* 2*fG(m,R0,0) * 1/(2*pi^2)
        % Afi = K * 2*fG(m,R0,0)
        % => K = pi*mi0*I / (2*pi^2) = mi0*I/(2*pi)
      endif
    endfor
    A0        = mi0*I/(2*pi) * A0;  % poračunamo vse konstante izven zanke
    A0        = A0 * cos(fi);       % vks2kks(A0)
  endif
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B0 = zanka(T, Ri,Zi, Rz, k)
  % magnetno polje v tocki T(r,z)
  % B0 = I * zanka(T, Ri,Zi, d, 1) !!!
  % tokovne zanke z radijem Ri na visini Zi
  global mi0
  [r fi z]  = kks2vks(T);
  [R0 m] = R0_m(r,z,Ri,Zi);
  while (m == 1)                           % tocka v osi zice?
    printf("ZANKA: m = %d zanke (%1.2e, %1.2e) v tocki (%1.2e, %1.2e)\n", m, Ri, Zi, r, z);
    r=r+Rz;                             % racunaj kot tocko na robu zice
    [R0 m] = R0_m(r,z,Ri,Zi);           % nov R0, m
  endwhile

  Br = 0; Bz = 0;
  % tocka T v VKS;  tocka na zanki Tz v VKS;  k .. mnozenje s konstanto
  if (r > 0)
    [K,E] = ellipke(m);               % Eliptic integrala
    Br = (Zi-z)/(r*R0) * (2*K - (2-m) * E/(1-m));
    Bz = (1/R0)*(2*K - (2 - m*(r+Ri)/r) * E/(1-m));
    B0 = [Br*cos(fi); 0; Bz];           % [Bx; 0; Bz]
    if(k) B0 = mi0/(4 *pi) * B0; endif  % k .. poracunam konstante?
  else
    Bz = Ri^2 / ((Ri^2+(z-Zi)^2)^(3/2)); % Br = 0
    B0 = [0; 0; Bz];
    if(k) B0 = mi0/2 * B0; endif        % k .. poracunam konstante?
  endif
endfunction

% ----------------------------------------------------------------------
function B0 = B0L(T, data)
%/**
% * Gostota magnetnega pretoka v tocki (x,y,z)
% * ELMG str. 201
% */
  global mi0
  [I, l, Rout, Rin, Zoff, d] = vec2var(data);
  [RR ZZ NN] = osiZank(data);

  B0  = zeros(3,1);
  for n = 1:NN;                   % po vseh zankah
  Ri = RR(n);   Zi = ZZ(n);
      B0 += zanka(T, Ri,Zi, d, 1);  % Izracun B0 v valjnem KS (r,z)
  endfor
  B0 = I*B0;                        % tok v vsaki zanki
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Izracun enojnega sloja
function KMM = calcKMM(dataMM, dataL)
  global mi0
  global ex ey ez
  [mi, lMM, rMM, Zoff, Ns, dl, NsO, dlO] = vec2var(dataMM);
  NsMM              = Ns + 2*NsO;
  efi               = ey;                   % za nase tocke!

  alphaK            = zeros(NsMM, NsMM);
  a                 = zeros(NsMM, 1);
  KMM               = zeros(NsMM, 1);

  for jj = 1:NsMM
% ----- izracun desne strani sistema enacb -----------------------------
    [Tj nj dlj]     = TnMM(jj);             % tocka in normala j-te tocke
    B0              = B0L(Tj, dataL);       % gostota mag. polja tuljave na robu MM
    tj              = cross(nj,efi);        % tangenta površine v točki Tj
    bj              = beta(mi);             % razmerje beta
    B0t             = dot(B0,tj);           % B0t(Tj) = B0(Tj) * tj (skalarno)
    a(jj)           = 2*bj/mi0 * B0t;

% ----- izracun leve strani sistema enacb ------------------------------
    for ii = 1:NsMM % seštevamo vplive vseh i-tih zank (tokov) na točko Tj
      [Ti ni dli] = TnMM(ii);               % tocka i-te zanke
      [ri fi zi]  = kks2vks(Ti);            % tocka na robu v valjni KS, za dimenzije zanke
      if ii==jj
        B1          = [0; 0; mi0/(4*pi*ri)*log(16*ri/dli)]; %Bz; Br = 0;
      else
        B1          = zanka(Tj, ri,zi, dlj, 1);  % i-ti prispevek B1
      endif
      B1t           = dot(B1, tj);          % TANGENTA v tocki Tj !!!
      B1t_K         = dlj * B1t;            % i-ti tok Ii = Ki*dli; 
      alphaK(jj,ii) = 2*bj/mi0 * B1t_K;
      if ii==jj
        alphaK(jj,ii) += 1;                 % B0
      endif
    endfor  % po vseh i-jih
  endfor  % po vseh j-jih
  KMM = linsolve (alphaK, a);
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function A1 = A1MM(T, KMM)
%/**
% * Gostota magnetnega pretoka v tocki (x,y,z)
% * ELMG str. 201
% */
  global mi0 dataMM
  % @param T an array argument; T = [x y z];
  % 	točka v kateri računamo polje
  % @param dim an array argument;
  %   data = [mi, l, Ri, Zoff, Ns, dl, NsO, dlO];
  % 	mi: ...... permeabilnost MM
  % 	l ........ dolzina MM
  % 	Ri: ...... polmer MM
  %   Zoff ..... pozicija tuljave po Z (koordinati)
  % 	Ns ....... izracunano stevilo segmentov
  %   dl ....... izračunana dolžina segmenta
  % 	NsO ...... izracunano stevilo segmentov os. ploskve
  %   dlO ...... izračunana dolžina segmenta O
  [mi, l, Ri, Zoff, Ns, dl, NsO, dlO] = vec2var(dataMM);
  NsMM            = Ns + 2*NsO;
  A1              = 0;
  [r fi z]        = kks2vks(T);
  for ii = 1:NsMM
    [Ti ni dl]    = TnMM(ii);         % tocka, normala in dolzina segmenta
    [Ri f Zi]     = kks2vks(Ti);
    [R0 m]        = R0_m(r,z, Ri,Zi);
    % Afi = 2*pi*mi0*Yo * fG(m,R0,1)  % racunamo enacbo na str.: 201 -----
    %                    I = KMM*dl
    if (r > 0) && (m != 1)            % Potencial na osi
      A1          = A1 + KMM(ii)*dl * fG(m,R0,0);  % superpozicija ovojev, fG brez konstant
    endif
  endfor
  A1              = A1 * mi0/pi;    % poračunamo vse konstante izven zanke
  A1              = A1 * cos(fi);   % nazaj v KKS
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B1 = B1MM(T, KMM)
%/**
% * Gostota magnetnega pretoka v tocki (x,y,z)
% * ELMG str. 201 ?
% */
  global mi0 dataMM
  % @param T an array argument; T = [x y z];
  % 	točka v kateri računamo polje
  % @param dim an array argument;
  %   data = [mi, l, Ri, Zoff, Ns, dl, NsO, dlO];
  % 	mi: ...... permeabilnost MM
  % 	l ........ dolzina MM
  % 	Ri: ...... polmer MM
  %   Zoff ..... pozicija tuljave po Z (koordinati)
  % 	Ns ....... izracunano stevilo segmentov
  %   dl ....... izračunana dolžina segmenta
  % 	NsO ...... izracunano stevilo segmentov os. ploskve
  %   dlO ...... izračunana dolžina segmenta O
  [mi, l, Ri, Zoff, Ns, dl, NsO, dlO] = vec2var(dataMM);
  NsMM            = Ns + 2*NsO;
  [r f z]         = kks2vks(T);
  B1 = zeros(3,1);
  for ii = 1:NsMM
    [Ti ni dl]    = TnMM(ii);       % tocka, normala in dolzina segmenta
    [Ri fi Zi]    = kks2vks(Ti);
    B1 += zanka(T, Ri,Zi, dl, 1);   % Br, Bz
  endfor
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = WmL(data, KMM = 0)
  global mi0 dataMM dataL1 dataL2
  [I, l, Rout, Rin, Zoff, d] = vec2var(data);
  [RR ZZ NN] = osiZank(data);
  
  if(KMM)                          % z ali brez vpliva MM
    [mi, lMM, rMM, Zoff, Ns, dl, NsO, dlO] = vec2var(dataMM);
    NsMM  = Ns + 2*NsO;
  endif   
  
  W       = 0;
  %% tocke izracuna
  Rz = d/2;
  for n = 1:NN;                   % po vseh zankah
    Ri  = RR(n);
    Zi  = ZZ(n);
    T   = [Ri, 0, Zi];
  
    %% izvori - skupni A!!
    A   = 0;
    A  += A0L(T, dataL1);           % A tuljave i, brez MM
    A  += A0L(T, dataL2);           % A tuljave i, brez MM
    if(KMM)
      A  += A1MM(T, KMM);             % A tuljave i z vplivom MM
    endif

    %% integral J*A dv
    % J  = I/(pi*Rz^2);                % gostota toka
    % dv = pi*Rz^2 dl
    % int dl = li
    li = 2*pi*Ri;                     % srednja dolzina
    W += I*li*A;
  endfor
endfunction

function W = getWm(KMM = 0)
  global dataMM dataL1 dataL2
  W   = WmL(dataL1, KMM);
  W  += WmL(dataL2, KMM);
endfunction

% ----------------------------------------------------------------------
function L = lastnaL(data, KMMi = 0)
%/**
% * Lastna in medsebojna induktivnost tuljave
% * ELMG str. 240
% */
  global mi0 dataMM
  [I, l, Rout, Rin, Zoff, d] = vec2var(data);
  [RR ZZ NN] = osiZank(data);
  Rz = d/2;
  
  if(KMMi)                          % z ali brez vpliva MM
    [mi, lMM, rMM, Zoff, Ns, dl, NsO, dlO] = vec2var(dataMM);
    NsMM  = Ns + 2*NsO;
  endif   
  
  L     = 0;
  %% tocke izracuna
  for n = 1:NN;                     % po vseh zankah
    Ri = RR(n);   Zi = ZZ(n);
    Ti = [Ri, 0, Zi];
  
    %% izvori
    Ai  = A0L(Ti, data);            % A tuljave i, brez MM
    if(KMMi)
      Ai += A1MM(Ti, KMMi);         % A tuljave i z vplivom MM
    endif

     %% integral J*A dv
    % J  = I/(pi*Rz^2);             % gostota toka
    % dv = pi*Rz^2 dl
    % int dl = li
    li = 2*pi*Ri;                   % srednja dolzina
    L += li*Ai;
  endfor
  L = 1/I * L;                      % L = 1/I^2 * int J*A dv 
endfunction

% **********************************************************************
% Pomozne funkcije
% **********************************************************************
function [R Z N] = osiZank(data)
  [I, l, Rout, Rin, Zoff, d] = vec2var(data);
  Rsum   = (Rin+d/2):d:Rout-d/2;
  Zsum   = (Zoff - l/2 + d/2):d:(Zoff+ l/2 - d/2);
  [R Z] = meshgrid(Rsum, Zsum);
  N = numel(R);
endfunction

function [T n dl] = TnMM(i)
global dataMM
% TnMM vrne koordinate sredisca in normalo i-tega segmenta na plascu
% funkcija TnMM nam v celoti opisuje geometrijo mehko magnetika
% podatke o dimenzijah pridobi iz globalnih spremenljivk
  [mi, l, R, Zoff, N, dl, NO, dlO] = vec2var(dataMM);
  Nsum      = (N+2*NO);
  z0        = l/2;
  M         = [0; 0; Zoff];                   % vektor zamika
  if(Nsum != 1)
    if (i <= NO)
      T     = [(i-.5)*dlO; 0; z0];
      n     = [0; 0; 1];
      dl    = dlO;
    elseif (i > N + NO)
      T     = [(R -(i-(N + NO)-.5)*dlO); 0; -z0];
      n     = [0; 0; -1];
      dl    = dlO;
    else
      T     = [R; 0; z0-(i-NO-.5)*dl];
      n     = [1; 0; 0];
    endif
  else
    T       = [R; 0; l/4];   n = [1; 0; 0];
  endif  
  T         = T + M;
endfunction

function b = beta(mi)
  global mi0
  b         = (mi - mi0)/(mi + mi0);
endfunction

% ----- function -------------------------------------------------------
function [s t T] = timestring(tp)
    t = now;
    if nargin == 0
        tp = t;
    end
    c = datevec(t);
    s = datestr (c, 0);

    dif = t-tp;
    cdif = datevec(dif);
    T = cdif(4:6);
endfunction


% ----- function -------------------------------------------------------
% Vector to variables
function [var1 var2 var3 var4 var5 var6 var7 var8 var9]  = vec2var(vec)
  vec(size(vec,2)+1:9) = 0;
  var1 = vec(1);
  var2 = vec(2);
  var3 = vec(3);
  var4 = vec(4);
  var5 = vec(5);
  var6 = vec(6);
  var7 = vec(7);
  var8 = vec(8);
  var9 = vec(9);
endfunction

% Matrix to vector
function vec = mat2vec(M)
  vec(1) = M(1);  vec(2) = M(4);  vec(3) = M(7);
  vec(4) = M(2);  vec(5) = M(5);  vec(6) = M(8);
  vec(7) = M(3);  vec(8) = M(6);  vec(9) = M(9);
endfunction

% Matrix to vector
function M = vec2mat(vec)
  M(1,1) = vec(1);  M(1,2) = vec(4);  M(1,3) = vec(7);
  M(2,1) = vec(2);  M(2,2) = vec(5);  M(2,3) = vec(8);
  M(3,1) = vec(3);  M(3,2) = vec(6);  M(3,3) = vec(9);
endfunction

% ----- function -------------------------------------------------------
% (x, y, z) to (rho, phi, z)
function [r, fi, z]  = kks2vks(T)
  r = sqrt(T(1)^2+T(2)^2);
  fi = atan2(T(2), T(1));
  z = T(3);
endfunction

% ----- function -------------------------------------------------------
% (rho, phi, z) to (x, y, z)
function [x, y, z]  = vks2kks(T)
  x = T(1)*sin(T(2));
  y = T(1)*cos(T(2));
  z = T(3);
endfunction


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                      %
% GRAFIČNE FUNKCIJE                                                    %
%                                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function narisiTuljavo(data, ovoji=0)
  [I, l, Rout, Rin, Zoff, d] = vec2var(data);
  Z0 = l/2;
  M = [0; Zoff];
  T0A = [-Rin;  -Z0] + M;
  T0B = [ Rin;  -Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');
  T0A = [-Rin;   Z0] + M;
  T0B = [ Rin;   Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');

  T0A = [-Rin;  -Z0] + M;
  T0B = [-Rout; -Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');
  T0A = [ Rin;  -Z0] + M;
  T0B = [ Rout; -Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');

  T0A = [-Rin;   Z0] + M;
  T0B = [-Rout;  Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');
  T0A = [ Rin;   Z0] + M;
  T0B = [ Rout;  Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');

  T0A = [-Rin;  -Z0] + M;
  T0B = [-Rin;   Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');
  T0A = [-Rout; -Z0] + M;
  T0B = [-Rout;  Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');

  T0A = [ Rin;  -Z0] + M;
  T0B = [ Rin;   Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');
  T0A = [ Rout; -Z0] + M;
  T0B = [ Rout;  Z0] + M;
  line([T0A(1) T0B(1)],[T0A(2) T0B(2)],'LineWidth',1,'Color','black');
  
  if(ovoji)
    [Rsum Zsum] = osiZank(data);
    [Rosi Zosi] = meshgrid(Rsum, Zsum);
    plot(Rosi, Zosi, "or");
  endif
endfunction

function narisiMM(data)
  Z0 = data(2)/2;
  M = [0; data(4)];
  rMM  = data(3);
  TA = [-rMM; -Z0] + M;
  TB = [ rMM; -Z0] + M;
  line([TA(1) TB(1)],[TA(2) TB(2)],'LineWidth',2,'Color','black');
  TA = [-rMM;  Z0] + M;
  TB = [ rMM;  Z0] + M;
  line([TA(1) TB(1)],[TA(2) TB(2)],'LineWidth',2,'Color','black');
  TA = [-rMM; -Z0] + M;
  TB = [-rMM;  Z0] + M;
  line([TA(1) TB(1)],[TA(2) TB(2)],'LineWidth',2,'Color','black');
  TA = [ rMM; -Z0] + M;
  TB = [ rMM;  Z0] + M;
  line([TA(1) TB(1)],[TA(2) TB(2)],'LineWidth',2,'Color','black');
endfunction

function narisiWm(data, Wm)
[Rsum Zsum] = osiZank(data)
endfunction
