function Gauss()
    Aufgabenteil = 'a'
    aufgabenteil_a();
    
    Aufgabenteil = 'b'
    %aufgabenteil_b();
    
    Aufgabenteil = 'c'
    %aufgabenteil_c();
end

function aufgabenteil_a()
    % A und b wie in Aufgabe 1
    A = [2 3 -1; 3 4.5 -2; 4 2 0];
    b = [4; -2; 3];
    % Wende Gauss mit Spaltenpivotsuche an, da sonst ein Element == 0 wird
    [A, b] = doGauss(A,b,1);
    % LR-Zerlegung von A
    A
    % L?sungsvektor x des LGS Ax=b
    b
end

function aufgabenteil_b()
    for l=6:20
        A = [10^(-l) 1; 1 1];
        b = [1; 2];
        % L?se Ax=b ohne und mit Spaltenpivotsuche
        [C, ohne_pivot] = doGauss(A, b, 0);
        [E, mit_pivot] = doGauss(A, b, 1);
        % Genauigkeit der Eingabe
        l
        % L?sungsvektor x des LGS Ax=b ohne Spaltenpivotsuche
        ohne_pivot
        % L?sungsvektor x des LGS Ax=b mit Spaltenpivotsuche
        mit_pivot
    end
end

function aufgabenteil_c()
    % beschr?nke MATLAB auf einen physischen Kern
    maxNumCompThreads(1);
    % Gr??e der quadratischen Matrix A
    n = [5;10;50;100;1000];
    for i=1:5;
        n(i)
        % erzeuge A als nxn-Matrix mit pseudozuf?lligen Elementen
        A = rand( n(i) );
        % erzeuge b als Vektor mit n pseudozuf?lligen Elementen
        b = rand( n(i), 1 );
        tic;
        [A, b] = doGauss(A, b, 1);
        % Zeit, die unser Algorithmus f¸r die L?sung des LGS ben?tigt
        t_we = toc
        tic;
        A\b;
        % Zeit, die die MATLAB-Implementation ben?tigt
        t_matlab = toc
    end
end

% F¸hre eine LR-Zerlegung durch und berechne den L?sungsvektor
% @param    Matrix      A               Die Matrix des zu l?senden LGS
% @param    Vektor      b               Der L?sungsvektor des zu l?senden
%                                       LGS
% @param    boolean     pivot_flag      ==1: f¸hre Spaltenpivotsuche durch
%                                       (sonst nicht)
% @return   Matrix      A               LR-Zerlegung von A
% @return   Vektor      b               L?sungsvektor x des LGS Ax=b
function [A, b] = doGauss( A, b, pivot_flag )
    [A, P] = LU_decomposition(A, pivot_flag);
    b = compute(A, b, P);
end

% Hilfsfunktion f¸r die LR-Zerlegung einer Matrix
% @param    Matrix      A               Die Matrix, die LR-zerlegt werden soll
% @param    boolean     pivot_flag
% @return   Matrix      A               LR-Zerlegung von A
% @return   Matrix      P               Permutationsmatrix
function [A, P] = LU_decomposition(A, pivot_flag)
    [m, n] = size(A);
    % erzeuge eine neue Permutationsmatrix P = I
    P = diag(ones(1,n));
    % betrachte i-te Zeile
    for i=1:(n-1)
        % soll eine Spaltenpivotsuche duruchgef¸hrt werden?
        if (pivot_flag == 1)
            [A, P] = column_pivot(A, P, i, n);
        end
        % betrachte j-te Zeile
        for j=(i+1):n
            % nach Satz 3.3 m¸ssen alle A(i,i)~=0 sein
            if (A(i,i)~=0)
                l = A(j,i)/A(i,i);
                % betrachte k-te Spalte
                for k=(i+1):n
                    A(j,k) = A(j,k) - l*A(i,k);
                end
                % L-Matrix
                A(j,i) = l;
            % sonst wird ein Fehler zur¸ckgegeben
            % die Ausgaben durch die Hauptfunktion finden weiterhin statt
            else
                error = 'ACHTUNG ACHTUNG! A(';
                error = cat(2,error,num2str(i));
                error = cat(2,error,',');
                error = cat(2,error,num2str(i));
                error = cat(2,error,')=0, vertrauen Sie den folgenden Ergebnissen nicht!')
                break
            end
        end
    end
end

% F¸hre eine Spaltenpivotsuche durch
% @param    Matrix      A               Die zu bearbeitende Matrix
% @param    Matrix      P               Das Produkt der bisherigen
%                                       Permutationsmatrizen
% @param    Integer     i               Die zu bearbeitende Spalte
% @param    Integer     n               Gr??e der Matrix A
% @return   Matrix      A               Die spaltenpivotdurchsuchte Matrix
% @return   Matrix      P               P^(i)*P
function [A, P] = column_pivot(A, P, i, n)
    max_value = 0;
    max_pos   = 0;
    % betrachte j-te Zeile
    for j=i:n
        if (abs(A(j,i)) > max_value)
            % speicher den Wert und die Zeile des optimalen Tauschpartners
            max_value = abs(A(j,i));
            max_pos   = j;
        end
    end
    % vertausche i-te und max_pos-te Zeile in der Matrix A
    A([i,max_pos],:) = A([max_pos,i],:);
    % erzeuge eine neue Permutationsmatrix P^(i) = I
    P_i = diag(ones(1,n));
    % vertausche i-te und max_pos-te Zeile in der Matrix P^(i)
    P_i([i,max_pos],:) = P_i([max_pos,i],:);
    % transponiere P^(i) um an eine Spaltenvertauschung zu gelangen
    P_i = transpose(P_i);
    % berechne P^(i)*P
    P = P_i * P;
end

% Berechne den L?sungsvektor des LGS PAx=Pb
% @param    Matrix      A
% @param    Vektor      b
% @param    Matrix      P
% @return   Vektor      b               L?sungsvektor x des LGS PAx=Pb
function b = compute(A, b, P)
    % wende die Permutationsmatrix P auf den Vektor b an, speicher das
    % Resultat in b
    b = P*b;
    % l?se Ly=b durch vorw?rtseinsetzen, speicher y in b
    b = forward_insert(A,b);
    % l?se Rx=b durch r¸ckw?rtseinsetzen, speicher x in b
    b = backward_insert(A,b);
end

% L?se das LGS Ly=b durch vorw?rtseinsetzen
% @param    Matrix      A
% @param    Vektor      b
% @result   Vektor      b               L?sungsvektor y des LGS Ly=b
function b = forward_insert(A,b)
    [m, n] = size(A);
    % betrachte i-te Zeile
    for i=2:n
        % betrachte j-te Spalte
        for j=1:(i-1)
            b(i) = b(i) - (A(i,j)*b(j));
        end
    end
end

% L?se das LGS Rx=b durch r¸ckw?rtseinsetzen
% @param    Matrix      A
% @param    Vektor      b
% @result   Vektor      b               L?sungsvektor x des LGS Rx=b
function b = backward_insert(A,b)
    [m, n] = size(A);
    % betrachte (n-i)-te Zeile 
    for i=0:(n-1)
        % betrachte (n-j)-te Spalte
        for j=0:(i-1)
            b(n-i) = b(n-i) - (A(n-i, n-j)*b(n-j));
        end
        b(n-i) = b(n-i)/A(n-i,n-i);
    end
end