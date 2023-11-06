%Glowna program do sprawozdania
clear;
%zad 1.

%Parametry
MaxSize=500;


%{ 
KOMENTARZ MULTILINE
%}

%{
for i=1:1:MaxSize
            disp('iteration ');
            disp(i);
           
            %%%%% ODKOMENTOWAC W ZALEZNOSCI OD PODPUNKTU !!!%%%%
            %A = vander(1:i); % A) Vandermonde
            %A = pascal(i);   % B) Pascal
            %A = hilb(i);     % C) Hilbert    
            b = A*ones(size(A,1),1); % dla podpuktow A B C
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            if(det(A) == 0)
                display('RANK PROBLEM')
            end
            if(isIllConditioned(decomposition(A)))
                display('Badly condictioned') 
            end
            % x_ = wynik przyblizony, x = wynik dokladny
        
            %Uruchomienie funkcji i zliczenie czasu
            tic;
            x_= F_Gaus(A,b);
 
            t_arr(i,2)=toc; 
            t_arr(i,1)=i;
            
            %Wynik dokładny
            x = A\b;
            
            det_A_arr(i,1)=i;
            det_A_arr(i,2)=det(A);
            
            cond_A_arr(i,1)=i;
            cond_A_arr(i,2)=cond(A);
        
            error(i,1)=i;
            error(i,2) = sum(sqrt((x-x_).^2));
            
            error_r(i,1)=i; 
            error_r(i,2) = sum(sqrt(((A*x_)-b).^2));
     

end

ti= tiledlayout(2,3);
title(ti, 'Vandermonde');

nexttile;
plot(t_arr(:,1),t_arr(:,2));
xlabel('Number of variables');
ylabel('time');

nexttile;
plot(error(:,1),error(:,2));
xlabel('Number of variables');
ylabel('L2 norm error');

nexttile;
plot(error_r(:,1),error_r(:,2));
xlabel('Number of variables');
ylabel('Rest L2 norm error');

nexttile;
plot(det_A_arr(:,1),det_A_arr(:,2));
xlabel('Number of variables');
ylabel('Det A');

nexttile;
plot(cond_A_arr(:,1),cond_A_arr(:,2));
xlabel('Number of variables');
ylabel('Cond A');

%}




%%%%%%%%%%%%%%%%%%%%%%%%%% PODPUNKT D %%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:1:MaxSize
            disp('iteration ');
            disp(i);
            %Przygotowanie macierzy
        
            % PODPUNKT D
            b= rand(i,1);
            [Q,R] = qr(b);
            rand_triu_matrix=triu(rand(i));
            for j=1:1:10
                %disp('k ');
                %disp(j);
                A = Q*(rand_triu_matrix^j);
         
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
                if(det(A) == 0)
                    display('RANK PROBLEM')
                end
                if(isIllConditioned(decomposition(A)))
                    display('Badly condictioned') 
                end
                % x_ = wynik przyblizony, x = wynik dokladny
            
                %Uruchomienie funkcji i zliczenie czasu
                tic;
                x_= F_Gaus(A,b);
     
                t_arr(i,j+1)=toc; 
                t_arr(i,1)=i;
                
                %Wynik dokładny
                x = A\b;
                
                det_A_arr(i,1)=i;
                det_A_arr(i,j+1)=det(A);
                
                cond_A_arr(i,1)=i;
                cond_A_arr(i,j+1)=cond(A);
            
                error(i,1)=i;
                error(i,j+1) = sum(sqrt((x-x_).^2));
                
                error_r(i,1)=i; 
                error_r(i,j+1) = sum(sqrt(((A*x_)-b).^2));
            end
     

end

ti= tiledlayout(2,3);
title(ti, 'Vandermonde');

nexttile;
plot(t_arr(:,1),t_arr(:,2));
plot(t_arr(:,1),t_arr(:,3));
plot(t_arr(:,1),t_arr(:,4));
plot(t_arr(:,1),t_arr(:,5));
plot(t_arr(:,1),t_arr(:,6));
plot(t_arr(:,1),t_arr(:,7));
plot(t_arr(:,1),t_arr(:,8));
plot(t_arr(:,1),t_arr(:,9));
plot(t_arr(:,1),t_arr(:,10));
plot(t_arr(:,1),t_arr(:,11));
xlabel('Number of variables');
ylabel('time');

nexttile;
plot(error(:,1),error(:,2));
plot(error(:,1),error(:,3));
plot(error(:,1),error(:,4));
plot(error(:,1),error(:,5));
plot(error(:,1),error(:,6));
plot(error(:,1),error(:,7));
plot(error(:,1),error(:,8));
plot(error(:,1),error(:,9));
plot(error(:,1),error(:,10));
plot(error(:,1),error(:,11));
xlabel('Number of variables');
ylabel('L2 norm error');

nexttile;
plot(error_r(:,1),error_r(:,2));
plot(error_r(:,1),error_r(:,3));
plot(error_r(:,1),error_r(:,4));
plot(error_r(:,1),error_r(:,5));
plot(error_r(:,1),error_r(:,6));
plot(error_r(:,1),error_r(:,7));
plot(error_r(:,1),error_r(:,8));
plot(error_r(:,1),error_r(:,9));
plot(error_r(:,1),error_r(:,10));
plot(error_r(:,1),error_r(:,11));
xlabel('Number of variables');
ylabel('Rest L2 norm error');

nexttile;
plot(det_A_arr(:,1),det_A_arr(:,2));
plot(det_A_arr(:,1),det_A_arr(:,3));
plot(det_A_arr(:,1),det_A_arr(:,4));
plot(det_A_arr(:,1),det_A_arr(:,5));
plot(det_A_arr(:,1),det_A_arr(:,6));
plot(det_A_arr(:,1),det_A_arr(:,7));
plot(det_A_arr(:,1),det_A_arr(:,8));
plot(det_A_arr(:,1),det_A_arr(:,9));
plot(det_A_arr(:,1),det_A_arr(:,10));
plot(det_A_arr(:,1),det_A_arr(:,11));
xlabel('Number of variables');
ylabel('Det A');

nexttile;
plot(cond_A_arr(:,1),cond_A_arr(:,2));
plot(cond_A_arr(:,1),cond_A_arr(:,3));
plot(cond_A_arr(:,1),cond_A_arr(:,4));
plot(cond_A_arr(:,1),cond_A_arr(:,5));
plot(cond_A_arr(:,1),cond_A_arr(:,6));
plot(cond_A_arr(:,1),cond_A_arr(:,7));
plot(cond_A_arr(:,1),cond_A_arr(:,8));
plot(cond_A_arr(:,1),cond_A_arr(:,9));
plot(cond_A_arr(:,1),cond_A_arr(:,10));
plot(cond_A_arr(:,1),cond_A_arr(:,11));
xlabel('Number of variables');
ylabel('Cond A');



