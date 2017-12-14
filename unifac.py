import numpy as np

U  = 5

% Este es el número de moléculas en la mezcla
m = 2    

% Este es el número de grupos funcionales en la mezcla
%g = 3   
g = 7

%T = 331.15  % K
%T = 328
%     Etanol - n-Hexano
%xj = [0.332 , 0.668]
%xj = [0.383 , 0.617]

%--------------------------------------------------------
% Agua - Isoamil alcohol - ácido acético
%     H2O CH3 CH2 CH  OH  COOH  COOCH3
% v1 = [1   0   0   0   0   0     0]'; % Agua
% v2 = [0   2   2   1   0   0     1]'; % Isoamil acetato
% v3 = [0   1   0   0   0   1     0]'; % Ácido acético

% v = [v1' ; v2' ; v3']';
% v = [v1' ; v3']';
v
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Agua - Isoamil acetato - ácido acético
%     H2O     CH3    CH2    CH     OH    COOH   COOCH3
R = [0.9200 0.9011 0.6744 0.4469 1.0000 1.3013  1.9031]';
Q = [1.4000 0.8480 0.5400 0.2280 1.2000 1.2240  1.7280]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Agua - Isoamil alcohol - Ácido acético
%    H2O     CH3     CH2     CH      OH      COOH    COOCH3
a = [0       300     300     300    -229.1  -14.09   72.8700;...
     1318    0       0       0       986.5   663.5   232.100;...
     1318    0       0       0       986.5   663.5   232.100;...
     1318    0       0       0       986.5   663.5   232.100;...
     353.5   156.4   156.4   156.4   0       199     101.100;...
    -66.17   315.3   315.3   315.3  -151     0      -256.300;...
     200.800 114.800 114.800 114.800 245.400 660.200 0      ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = exp(-a./T);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : 1 : m

   r(:,j) = sum(R.*v(:,j));

end
r;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : 1 : m

   q(:,j) = sum(Q.*v(:,j));

end
q;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : 1 : m

   J(:,j) = r(1,j)*xj(1,j)/sum(r.*xj);

end
J;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j = 1 : 1 : m

   L(1,j) = q(1,j)*xj(1,j)/sum(q.*xj);

end
L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
li = 5.*(r - q) - (r - 1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lnYCi = log(J./xj) + 5.*q.*log(L./J) + li - (J./xj).*(sum(xj.*li));

%lnY1C = log(J(1,1)/xj(1,1)) + 5*q(1,1)*log(L(1,1)/J(1,1)) + li(1,1) - (J(1,1)/xj(1,1))*(xj(1,1)*li(1,1) + xj(1,2)*li(1,2))
%lnY2C = log(J(1,2)/xj(1,2)) + 5*q(1,2)*log(L(1,2)/J(1,2)) + li(1,2) - (J(1,2)/xj(1,2))*(xj(1,1)*li(1,1) + xj(1,2)*li(1,2))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeficiente de actividad residual del grupo (k)
% en la molecula (i) %%%%%%%%%%%%%%%%%%%%%%
% Fracción molar del grupo funcional (k)
% en la molecula (i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 1 : m      %Molécula (i)

   for k = 1 : 1 : g   %Grupo funcional (k)

      xg(k,i) = v(k,i)./sum(v(:,i));

   end

end

xg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 1 : m      %Molécula (i)

   for k = 1 : 1 : g   %Grupo funcional (k)

      Lg(k,i) = Q(k,1)*xg(k,i)/sum(Q.*xg(:,i));

   end

end

Lg;

%mor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 1 : m      %Molécula (i)

   for k = 1 : 1 : g   %Grupo funcional (k)

      ST(k,i) = sum(Lg(:,i).*A(:,k));

   end

end

ST = ST';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 1 : m      %Molécula (i)

   for k = 1 : 1 : g   %Grupo funcional (k)

      if i == 1
         STa(k,:) = (Lg(:,i)'.*A(k,:));
      elseif i == 2
         STa(k + g,:) = (Lg(:,i)'.*A(k,:));
      elseif i == 3
         STa(k + 2*g,:) = (Lg(:,i)'.*A(k,:));
      end

   end

end

STa;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1 : 1 : m      %Molécula (i)

   for k = 1 : 1 : g   %Grupo funcional (k)

      if i == 1
         lnTg(i,k) = Q(k,1).*(1 - log(ST(i,k)) - sum(STa(k,:)./ST(i,:)));
      elseif i == 2
         lnTg(i,k) = Q(k,1).*(1 - log(ST(i,k)) - sum(STa(k+g,:)./ST(i,:)));
      elseif i == 3
         lnTg(i,k) = Q(k,1).*(1 - log(ST(i,k)) - sum(STa(k+2*g,:)./ST(i,:)));
      end

   end

end
%lnT(1,k) = Q(k,1).*(1 - log(STg(1,k)) - sum(STga(k,:)./STg));

lnTg;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coeficiente de actividad residual del grupo (k)
% en la mezcla %%%%%%%%%%%%%%%%%%%%%%
% Fracción molar del grupo funcional (k)
% en la mezcla
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : 1 : m      %Molécula (i)

   STq(:,i) = sum(v(:,i)*xj(:,i));

end

STq = sum(STq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : 1 : g   %Grupo funcional (k)

   xs(k,:) = (sum(v(k,:).*xj))./(STq);

end

xs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1 : 1 : g   %Grupo funcional (k)

   Lgs(k,1) = Q(k,1)*xs(k,1)/sum(Q.*xs);

end

Lgs;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k = 1 : 1 : g   %Grupo funcional (k)

   STg(k,:) = sum(Lgs.*A(:,k));

end

STg = STg';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : 1 : g   %Grupo funcional (k)

   STga(k,:) = (Lgs'.*A(k,:));

end

STga;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1 : 1 : g   %Grupo funcional (k)

  lnT(1,k) = Q(k,1).*(1 - log(STg(1,k)) - sum(STga(k,:)./STg));

end

lnT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Coeficiente de actividad Residual

for i = 1 : 1 : m      %Molécula (i)

   lnYRi(:,i) = sum(v(:,i).*(lnT' - lnTg(i,:)'));

lnYRi;

%Coeficiente de actividad total
lnYi = lnYCi + lnYRi
Yi = np.exp(lnYi)
