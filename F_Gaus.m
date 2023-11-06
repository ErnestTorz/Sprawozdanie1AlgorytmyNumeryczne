%Ernest Torz 7853 informatyka niestacjo mgr 1 sem
%Współczynniki macierzy
function out = F_Gaus(WS,WW)
%WS Wyrazy wolne
%WW Wspolczynniki
    A = [WS WW]; %Budowanie macierzy
    n= size(A,1); %Ilość wierszy                                                                                    
    
    %Iteracja po kolumnch i wierszach
    
    %i - iteracja po przekatnej (n-1, poniewaz dla ostatniego elementu przekatnej dzialania sa niepotrzebne 
     for i=1:n-1
        %[M,MI]=max(abs(A(i:n,i))) %Znajdowanie max abs w danej kolumnie poniżej diagonali
        %MI=MI+i-1 %Kalkulacja indeksu wiersza
        %A([i MI],:)=A([MI i],:) %Zamiana wierszy
        for j=i+1:n %j - iteracja po elementach w dół, od przekatnej
            m = A(j,i)/A(i,i); %Wzór na m dla danej komurki macierzy
            A(j,:) = A(j,:) - m*A(i,:); % Odejmowanie wierszy
        end
    end
    
    x = zeros(n,1); %Macierz dla rozwiazan
    x(n) = A(n,n+1)/A(n,n); %Wyliczamy zmienna dla ostatniego wiersza 
     for i=n-1:-1:1 %Iteracja od przedostatniego wiersza do pierwszego
        summ = 0;
        for j=i+1:n 
            summ = summ + A(i,j)*x(j) ;
        end
        x(i) = (A(i,n+1) - summ)/A(i,i); %Od wyrazu wolnego odejmujemy sume i dzielmi przez element przekatnej 
     end
     out = x;
end
