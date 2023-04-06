# studio-di-un-profilo-termico-con-matlab

Richieste dell’esercitazione:
Studiare la soluzione numerica di un problema di diffusione del calore in una sbarretta di lunghezza 1 metro e
spessore trascurabile i cui estremi vengono posti alle temperature di 0oC e di 100oC. Il profilo di temperatura
iniziale è dato da T(x,0) = 100*sin(πx/L).
In particolare si dovranno sviluppare i seguenti punti:
1) Analisi di convergenza della soluzione numerica alla soluzione analitica al variare del numero di intervalli
di discretizzazione spaziale N=(10, 20, 30)
2) Analisi di convergenza della soluzione numerica alla soluzione analitica al variare del numero di intervalli
di discretizzazione spaziale durante il transitorio N=(10, 20, 30)
3) Analisi della soluzione numerica al variare della lunghezza della sbarretta L=(0.5, 1, 1.5 m)
4) Analisi della soluzione numerica al variare del gradiente termico T = (10, 100, 1000oC)
5) Analisi della soluzione numerica al variare della diffusività termica del materiale (Argento, Ferro, Vetro)

codice Matlab:
clear
L = input('inserire la lunghezza della sbaretta (0.5, 1, 1.5) [m] = ')
num_nodi = input('inserire il numero di intervalli di discretizzazione spaziale (10, 20 o 30) = ') + 1
v_nodi = linspace(0, L, num_nodi)    %x
t_in(1:num_nodi) = 100*sin( pi * v_nodi / L);  % T iniziali 
t_in(end) = input('inserire la temperatura dell estremo destro (10, 100, 1000) [°C] = ') %la T del nodo finale è fissa
      a = 1.656 * 10^(-4); f = 1.88 * 10^(-5); v = 3.4 * 10^(-7);  %diffusività termica di argento, ferro e vetro
      materiale = input('scegliere il materiale tra argento, ferro e vetro. \n Selezionare tra: \n a = argento \n f = ferro \n v = vetro \n')
      if materiale == a
          D = a
      elseif materiale == f
          D = f
      elseif materiale == v
          D = v
      end
      
%soluzione ANALITICA
t_an = ( t_in(end) - t_in(1) ) * v_nodi / L

%grafico soluzione analitica
figure, plot(v_nodi, t_an, 'Marker', 'd', 'Color', 'r') 
title('T nei nodi alla fine, soluzione analitica'), xlabel('nodi'), ylabel('temperature finali');


%grafico iniziale delle T
figure, plot(v_nodi, t_in, 'Marker', 'd', 'Color', 'r') 
title('T nei nodi all inizio'), xlabel('nodi'), ylabel('temperature iniziali');



%problema STAZIONARIO
t_new = t_in
figure('Color', 'w', 'WindowState', 'maximize')
subplot(2,2,1)
h = animatedline('Marker', 'd', 'Color', 'r'), title('T nei nodi a regime stazionario'),
xlabel('nodi'), ylabel('temperature');
eps = 0.001   %errore massimo concesso
err = eps + 1  %inizializzazione dei valori
n_iter = 0;    %n di iterazioni
 while (err > eps && n_iter < 700)  %diamo un limite massimo di iterazioni
    for i = 2 : (num_nodi - 1)      %si sposta lungo i nodi, senza tocccare le T degli estremi
        t_new(i) = (t_new(i-1) + t_new(i+1)) / 2;
        i = i + 1;
    end
    err = norm( t_new - t_an, inf );
    t_new;
    addpoints(h, v_nodi, t_new)
    drawnow
    n_iter = n_iter + 1;
    vettore_err(n_iter) = err;
    vettore_iter(n_iter) = n_iter;   %vettore delle iterazioni
    clearpoints(h)
    drawnow limitrate
 end
addpoints(h, v_nodi, t_new)
t_new

%problema TRANSITORIO
tempo_convergenza_tr = 0 % lo inizializzo a zero
delta_x = L / (num_nodi - 1)   % passo = Lunghezza / intervalli di discretizzazione(= num_nodi-1)
tau = (delta_x^ 2)/(3 * D)
M_inv = D * tau / (delta_x)^ 2
t_new_tr = t_in    % vettore delle T del transitorio
subplot(2,2,2)
h_tr = animatedline('Marker', 'o', 'Color', 'g'), title('T nei nodi durante il transitorio'),
xlabel('nodi'), ylabel('temperature');
eps_tr = 0.01   %errore massimo concesso
err_tr = eps_tr + 1  %inizializzazione dei valori
n_iter_tr = 0;   %n di iterazioni
 while (err_tr > eps_tr && n_iter_tr < 700)   %diamo un limite massimo di iterazioni
     for i = 2 : (num_nodi - 1)
        t_new_tr(i) = t_new_tr(i-1)*M_inv + t_new_tr(i+1)*M_inv + t_new_tr(i)*(1 - 2*M_inv) 
        i = i + 1;
     end
    err_tr = norm( t_new_tr - t_an, inf );
    t_new_tr;
    addpoints(h_tr, v_nodi, t_new_tr)
    drawnow
    n_iter_tr = n_iter_tr + 1;
    vettore_err_tr(n_iter_tr) = err_tr;
    vettore_iter_tr(n_iter_tr) = n_iter_tr;   %vettore delle iterazioni
    clearpoints(h_tr)
    drawnow limitrate
    tempo_convergenza_tr = tempo_convergenza_tr + tau 
end
addpoints(h_tr, v_nodi, t_new_tr)
t_new_tr

% grafico degli errori 
subplot(2, 2, 3), plot(vettore_iter, vettore_err, 'LineWidth',2), 
title('errore a regime stazionario'), xlabel('iterazioni'), ylabel('errore')
% grafico degli errori          
subplot(2, 2, 4), plot(vettore_iter_tr, vettore_err_tr, 'LineWidth',2),
title('errore durante il transitorio'), xlabel('iterazioni'), ylabel('errore')

display ('la diffusivita termica è:')
display (D)

display ('l intervallo di tempo è:')
display (tau)

display ('la durata del transitorio è :')
display (tempo_convergenza_tr)

