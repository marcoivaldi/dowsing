%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**    written by Lorenzo Giacometti, research center SUISM - Torino    **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%
%**                                                                     **%
%**                                                                     **%
%** Media, SD, mediana, ARV, varianza delle varie epoche delle varie    **%
%** azioni, dei vari file .sig importati e stampa stampa tabelle        **%
%** su file .csv (analisi in base alle varie bande di frequenza)        **%
%** tenendo fisso 1 soggetto alla volta ponendo in X i canali e         **%
%** in Y i task, tenendo fisso un task alla alla volta ponendo          **%
%** in X i canali e in Y i soggetti, tenendo fisso sempre il task       **%
%** ma ponendo in X i soggetti e in Y i canali suddividendo l'emisfero  **%
%** destro e quello sinistro.											**%
%** Il programma da l'opportunità di costruire più tabelle, una ogni 	**%
%** banda di frequenza inserita.										**%
%**                                                                     **%
%**                                                                     **%
%*************************************************************************%
%*************************************************************************%

clear;clc;

versione = 'ANALISI DATI Versione 4.1.3d';

%definizioni per l'inserimento dei valori delle bande
f_def = {'prima','seconda','terza','quarta','quinta','sesta','settima','ottava','nona','decima','undicesima','dodicesima','tredicesima','quattordicesima','quindicesima','sedicesima','diciassettesima','diciottesima','diciannovesima','ventesima'};

answer = 0;

%input n°file da importare
answer = inputdlg('Inserisci il numero di file da elaborare','N° file da elaborare',1);
answer = str2double(answer);
n_file = answer(1,1);       %n° di file impostato da inserire

answer = 0;

%input n°bande da analizzare
answer = inputdlg('Inserisci il numero di bande da analizzare','N° bande da analizzare',1);
answer = str2double(answer);
n_band = answer(1,1);       %n° di bande da inserire

%frequenze iniziali bande
fr_inf = inputdlg(f_def(1:n_band),'Inserimento valori inferiori delle bande',1);
fr_inf = str2double(fr_inf);

%frequenze finali bande
fr_sup = inputdlg(f_def(1:n_band),'Inserimento valori superiori delle bande',1);
fr_sup = str2double(fr_sup);

canali = 16;                %numero di canali da analizzare dal file sig
freq_camp = 512;            %frequenza di campionamento segnale 
max_dim = 55;               %impostazione della massima lunghezza della sezione in secondi
dim_sec = [25,25,25,55,55,25,25,55,55];               %impostazione delle lunghezze delle sezioni in secondi

n_sez = 9;                  %n° delle sezioni da analizzare

tasks = {'OA','OC','OCSA','BIPOC','BIPOA','MONOC','MONOA','CUBDIS','CUBORD'};

%ordine dei canali per emisfero, per l'analisi soggetto_x_canale 
ordine_DX = [5, 12, 6, 7, 13, 8, 14, 16];
ordine_SX = [1, 9, 2, 3, 10, 4, 11, 15];

%preallocazioni:
ris_var(n_band,n_file,n_sez,canali) = double(0);   
ris_arv(n_band,n_file,n_sez,canali) = double(0);   
ris_med(n_band,n_file,n_sez,canali) = double(0);   
ris_mean(n_band,n_file,n_sez,canali) = double(0);   
ris_devstd(n_band,n_file,n_sez,canali) = double(0);  
ris_mean_round(n_band,n_file,n_sez,canali) = double(0);  
ris_density(n_band,n_file,n_sez,canali) = double(0);  
ris_max_density(n_band,n_file,n_sez,canali) = double(0);  
ris_X_XII_density(n_band,n_file,n_sez,canali) = double(0);  

varianza(n_band,n_file,n_sez,canali,max_dim) = double(0);
arv(n_band,n_file,n_sez,canali,max_dim) = double(0);
fmean(n_band,n_file,n_sez,canali,max_dim) = double(0);
fmed(n_band,n_file,n_sez,canali,max_dim) = double(0);
devstd(n_band,n_file,n_sez,canali,max_dim) = double(0);


%etichette input sezioni
s_iniz = {'Occhi Aperti:','Occhi Chiusi:','Occhi Chiusi e Sguardo Alto:','Appoggio bipodalico ad occhi chiusi','Appoggio bipodalico ad occhi aperti','Appoggio monopodalico ad occhi chiusi','Appoggio monopodalico ad occhi aperti','Cubo Rubik disordine', 'Cubo Rubik ordine'};
s_fine = {'Occhi Aperti:','Occhi Chiusi:','Occhi Chiusi e Sguardo Alto:','Appoggio bipodalico ad occhi chiusi','Appoggio bipodalico ad occhi aperti','Appoggio monopodalico ad occhi chiusi','Appoggio monopodalico ad occhi aperti','Cubo Rubik disordine', 'Cubo Rubik ordine'};
def_in = {'18', '63', '108','000','000','000','000','153','228'};
def_fin = {'43', '88', '133','000','000','000','000','208','283'};

for f=1:1:n_file               %ciclo che scandisce i file

    PathLoad = 0;
    LoadIndex = 0;

    [InFile,PathLoad, LoadIndex] = uigetfile('*.sig','Caricamento file .sig'); 

    if LoadIndex == 0 %se path == 0 significa che il file non è stato caricato 
        break; % e termino in programma
    end

    files(f,:) = InFile;       %sistemare l'arrey di stringhe dei nomi dei file per la stampa

    h = fopen([PathLoad InFile],'r');
    sig = fread(h,[17 inf],'short');

    sig = sig';

    tmp1 = exist([PathLoad InFile(1:20) '.mat'],'file');        %controllo l'esistenza di un file con i range già salvati
    
    if tmp1 == 2
        
        % Construct a questdlg with three options
        choice1 = questdlg('è stato trovato un file dati relativo,vuoi caricarlo?','Caricamento file dati','Carica','Inserisci da zero','Controlla i dati','Inserisci da zero');
        
        % Handle response
        switch choice1
            case 'Carica'
                load([PathLoad InFile(1:20) '.mat'],'sec0','secf');
                %break			
            case 'Controlla i dati'
                load([PathLoad InFile(1:20) '.mat'],'sec0','secf');
               
                % Construct a questdlg with three options
               
                str_sec0 = num2str(sec0');
                str_secf = num2str(secf');
                
                asd = {['Secondi iniziali: ' str_sec0] , [' Secondi finali: ' str_secf]};
                choice2 = questdlg(asd,'Verifica file dati','Carica','Inserisci da zero','Inserisci da zero');
                
                % Handle response
                switch choice2
                    case 'Carica'                     
                        %break			
                    case 'Inserisci da zero'
                       
                        sec0 = inputdlg(s_iniz,'Inserimento valori iniziali delle sezioni',1,def_in);
                        sec0 = str2double(sec0);

                        %tempo fine in secondi %tempo inizio in secondi
                     
                        secf = inputdlg(s_fine,'Inserimento valori finali delle sezioni',1,def_fin);
                        secf = str2double(secf);

                        T0 = (sec0* freq_camp) +1;              %tempo inizio in campioni
                        Tf = (secf * freq_camp) +1;             %tempo fine in campioni

                        % n° di secondi da analizzare nelle varie fasi da confrontare
                        %dim = [(secf(1)-sec0(1)), (secf(2)-sec0(2)), (secf(3)-sec0(3))];    
                        dim = secf-sec0;
                        dim = dim';
                        
                        %salva file secondi
                        save([PathLoad InFile(1:20) '.mat'],'sec0','secf');

                       %break                   
                end
                      
               %break
            case 'Inserisci da zero'
             
                sec0 = inputdlg(s_iniz,'Inserimento valori iniziali delle sezioni',1,def_in);  %SISTEMARE...
                sec0 = str2double(sec0);

                %tempo fine in secondi %tempo inizio in secondi
                secf = inputdlg(s_fine,'Inserimento valori finali delle sezioni',1,def_fin);
                secf = str2double(secf);

                T0 = (sec0* freq_camp) +1;              %tempo inizio in campioni
                Tf = (secf * freq_camp) +1;             %tempo fine in campioni

                % n° di secondi da analizzare nelle varie fasi da confrontare
                %dim = [(secf(1)-sec0(1)), (secf(2)-sec0(2)), (secf(3)-sec0(3))];    
                dim = secf-sec0;
                dim = dim';
                
                %salva file secondi
                save([PathLoad InFile(1:20) '.mat'],'sec0','secf');

        end
        
       
    else
    
        sec0 = inputdlg(s_iniz,'Inserimento valori iniziali delle sezioni',1,def_in);
        sec0 = str2double(sec0);

        %tempo fine in secondi %tempo inizio in secondi
        
        secf = inputdlg(s_fine,'Inserimento valori finali delle sezioni',1,def_fin);
        secf = str2double(secf);

        T0 = (sec0* freq_camp) +1;              %tempo inizio in campioni
        Tf = (secf * freq_camp) +1;             %tempo fine in campioni

        % n° di secondi da analizzare nelle varie fasi da confrontare
        %dim = [(secf(1)-sec0(1)), (secf(2)-sec0(2)), (secf(3)-sec0(3))];    
        dim = secf-sec0;
        dim = dim';
        
        %salva file secondi
        save([PathLoad InFile(1:20) '.mat'],'sec0','secf');

    end
    
    for b=1:1:n_band
        
        %oltre al 5° grado l'uscita è sempre NaN
        
        [B,A] = butter(5,[fr_inf(b)/(freq_camp/2) fr_sup(b)/(freq_camp/2)]);    %imposta i dati per il filtro di butterworth
        sig_filt= filtfilt(B,A, sig);                       %filtraggio

        sig_filt = sig_filt';                               %trasposizione de sig filtrato per avere i canali come primo indice

        for s=1:1:n_sez                 %ciclo che scandisce le sezioni
            for c=1:1:canali            %ciclo che scandisce i canali

                clear P;
                for t=1:1:dim_sec(s)       	%ciclo che scandisce i secondi
                    %ARV:
                    sig_epoch= sig_filt(c,(sec0(s)+t-1)*freq_camp+1:(sec0(s)+t)*freq_camp);
                    arv(b,f,s,c,t)= sum(abs(sig_epoch))/length(sig_epoch);
                    
                    %calcolo frequenza media e frequenza mediana (MNF)
                    [fmean(b,f,s,c,t), fmed(b,f,s,c,t)]= fmed3cla(sig_epoch,freq_camp);

                    sig_epoch = sig_epoch - mean(sig_epoch);
                    
                    %PSD:
                    [P(t,:),fr] = psd(sig_epoch,freq_camp,freq_camp,boxcar(length(sig_epoch)),0);
                    
                    
                    %calcolo della varianza:
                    %varianza(b,f,s,c,t) = var(sig_filt(c,(sec0(s)+t-1)*freq_camp+1:(sec0(s)+t)*freq_camp));  %varianza di ogni epoca (1 secondo) del canale di una sezione
                    
                    %calcolo della deviazione standard
                    %devstd(b,f,s,c,t) = std(sig_epoch);

                end

                %ris_var(b,f,s,c) = mean(varianza(b,f,s,c,1:t));     %medio tutte le varianze della sezione s, del canale c, del file f
                ris_arv(b,f,s,c) = sum(arv(b,f,s,c,1:t));           %sommo tutti gli ARV della sezione s, del canale c, del file f
                ris_med(b,f,s,c) = mean(fmed(b,f,s,c,1:t));         %medio tutte le mediane della sezione s, del canale c, del file f
                ris_mean(b,f,s,c) = mean(fmean(b,f,s,c,1:t));       %medio tutte le medie della sezione s, del canale c, del file f
                %ris_devstd(b,f,s,c) = mean(devstd(b,f,s,c,1:t));    %medio tutte le deviazioni standard della sezione s, del canale c, del file f
                ris_mean_round(b,f,s,c) = round(ris_mean(b,f,s,c)); %arrotondo i valori di frequenza media in una nuova variabile con un solo decimale???.
                
                P = sum(P);
                P = P/dim_sec(s);
                ris_density(b,f,s,c) = P(ris_mean_round(b,f,s,c)+1);
                ris_max_density(b,f,s,c) = max(P);
                ris_X_XII_density(b,f,s,c) = P(9)+P(10)+P(11);
                
            end
        end
    end
end



files_original = files;                 %salvataggio originale dei nomi dei files
files = cellstr(files(:,3:6));          %estrazione nome soggetto dal nome del file
PathSave = uigetdir;                    %richiesta e salvataggio della path dove salvare gli output (matrici - file .csv)

%creazione delle cartelle per le diverse stampe
mkdir(PathSave,'CH_x_TSK');
mkdir(PathSave,'CH_x_SUBJ');
mkdir(PathSave,'SUBJ_x_CH_DX');
mkdir(PathSave,'SUBJ_x_CH_SX');
mkdir(PathSave,'SUBJ_x_TSK');

if n_band == 1    
    
    for b=1:1:n_band
        for f=1:1:n_file

            %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
            clear print_arv print_mean print_med

            %blocco l'incide della banda e quello del soggetto(file) per salvare task x canali 
            print_arv(1:n_sez,1:canali) = ris_arv(b,f,:,:);
            print_mean(1:n_sez,1:canali) = ris_mean(b,f,:,:);
            print_med(1:n_sez,1:canali) = ris_med(b,f,:,:);

            %traspongo le matrici per avere canali x task
            print_arv = print_arv';
            print_mean = print_mean';
            print_med = print_med';

            %stampa delle matrici su file .csv (separatore ,)
            csvwrite([PathSave '\CH_x_TSK\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_soggetto_' num2str(f) '.csv'],print_arv);
            csvwrite([PathSave '\CH_x_TSK\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_soggetto_' num2str(f) '.csv'],print_mean);
            csvwrite([PathSave '\CH_x_TSK\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_soggetto_' num2str(f) '.csv'],print_med);
        end       
    end

    for b=1:1:n_band
        for s=1:1:n_sez

            %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
            clear print_arv print_mean print_med

            %blocco l'incide della banda e quello del task per salvare canali x soggetti
            print_arv(1:n_file,1:canali) = ris_arv(b,:,s,:);
            print_mean(1:n_file,1:canali)= ris_mean(b,:,s,:);
            print_med(1:n_file,1:canali) = ris_med(b,:,s,:);

            %traspongo le matrici per avere canali x soggetti
            print_arv = print_arv';
            print_mean = print_mean';
            print_med = print_med';

            %stampa delle matrici su file .csv (separatore ,)
            csvwrite([PathSave '\CH_x_SUBJ\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_arv);
            csvwrite([PathSave '\CH_x_SUBJ\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_mean);
            csvwrite([PathSave '\CH_x_SUBJ\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_med);
        end       
    end

    for b=1:1:n_band
        for s=1:1:n_sez

            %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
            clear print_arv print_mean print_med print_arv_DX print_arv_SX

            %preallocazioni
            print_arv_DX(n_file,canali/2) = double(0);
            print_arv_SX(n_file,canali/2) = double(0);
            print_mean_DX(n_file,canali/2) = double(0);
            print_mean_SX(n_file,canali/2) = double(0);
            print_med_DX(n_file,canali/2) = double(0);
            print_med_SX(n_file,canali/2) = double(0);

            %blocco l'incide della banda e quello del task per salvare soggetti x canali
            print_arv(1:n_file,1:canali) = ris_arv(b,:,s,:);
            print_mean(1:n_file,1:canali) = ris_mean(b,:,s,:);
            print_med(1:n_file,1:canali) = ris_med(b,:,s,:);

            for i=1:1:8     %divisione in emisferi
                print_arv_DX(:,i) =  print_arv(:,ordine_DX(i));
                print_arv_SX(:,i) =  print_arv(:,ordine_SX(i));
                print_mean_DX(:,i) =  print_mean(:,ordine_DX(i));
                print_mean_SX(:,i) =  print_mean(:,ordine_SX(i));
                print_med_DX(:,i) =  print_med(:,ordine_DX(i));
                print_med_SX(:,i) =  print_med(:,ordine_SX(i));
            end

            %stampa delle matrici su file .csv (separatore ,)
            %emisfero DX
            csvwrite([PathSave '\SUBJ_x_CH_DX\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_arv_DX);
            csvwrite([PathSave '\SUBJ_x_CH_DX\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_mean_DX);
            csvwrite([PathSave '\SUBJ_x_CH_DX\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_med_DX);
            %emisfero SX
            csvwrite([PathSave '\SUBJ_x_CH_SX\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_arv_SX);
            csvwrite([PathSave '\SUBJ_x_CH_SX\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_mean_SX);
            csvwrite([PathSave '\SUBJ_x_CH_SX\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_med_SX);
        end       
    end

    for b=1:1:n_band
        for c=1:1:canali

            %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
            clear print_arv print_mean print_med

            %blocco l'incide della banda e quello del canale per salvare soggetti x task 
            print_arv(1:n_file,1:n_sez) = ris_arv(b,:,:,c);
            print_mean(1:n_file,1:n_sez) = ris_mean(b,:,:,c);
            print_med(1:n_file,1:n_sez) = ris_med(b,:,:,c);

            %stampa delle matrici su file .csv (separatore ,)
            csvwrite([PathSave '\SUBJ_x_TSK\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_arv);
            csvwrite([PathSave '\SUBJ_x_TSK\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_mean);
            csvwrite([PathSave '\SUBJ_x_TSK\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_med);
        end       
    end

else

    for b=1:1:n_band
    for f=1:1:n_file

        %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
        clear print_arv print_mean print_med

        %blocco l'incide della banda e quello del soggetto(file) per salvare task x canali 
        print_arv(1:n_sez,1:canali) = ris_arv(b,f,:,:);
        print_mean(1:n_sez,1:canali) = ris_mean(b,f,:,:);
        print_med(1:n_sez,1:canali) = ris_med(b,f,:,:);

        %traspongo le matrici per avere canali x task
        print_arv = print_arv';
        print_mean = print_mean';
        print_med = print_med';

        %stampa delle matrici su file .csv (separatore ,)
        csvwrite([PathSave '\CH_x_TSK\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_soggetto_' num2str(f) '.csv'],print_arv);
        csvwrite([PathSave '\CH_x_TSK\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_soggetto_' num2str(f) '.csv'],print_mean);
        csvwrite([PathSave '\CH_x_TSK\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_soggetto_' num2str(f) '.csv'],print_med);
    end       
    end

    for b=1:1:n_band
    for s=1:1:n_sez

        %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
        clear print_arv print_mean print_med

        %blocco l'incide della banda e quello del task per salvare canali x soggetti
        print_arv(1:n_file,1:canali) = ris_arv(b,:,s,:);
        print_mean(1:n_file,1:canali)= ris_mean(b,:,s,:);
        print_med(1:n_file,1:canali) = ris_med(b,:,s,:);

        %traspongo le matrici per avere canali x soggetti
        print_arv = print_arv';
        print_mean = print_mean';
        print_med = print_med';

        %stampa delle matrici su file .csv (separatore ,)
        csvwrite([PathSave '\CH_x_SUBJ\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_arv);
        csvwrite([PathSave '\CH_x_SUBJ\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_mean);
        csvwrite([PathSave '\CH_x_SUBJ\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_med);
    end       
    end

    for b=1:1:n_band
    for s=1:1:n_sez

        %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
        clear print_arv print_mean print_med print_arv_DX print_arv_SX

        %preallocazioni
        print_arv_DX(n_file,canali/2) = double(0);
        print_arv_SX(n_file,canali/2) = double(0);
        print_mean_DX(n_file,canali/2) = double(0);
        print_mean_SX(n_file,canali/2) = double(0);
        print_med_DX(n_file,canali/2) = double(0);
        print_med_SX(n_file,canali/2) = double(0);

        %blocco l'incide della banda e quello del task per salvare soggetti x canali
        print_arv(1:n_file,1:canali) = ris_arv(b,:,s,:);
        print_mean(1:n_file,1:canali) = ris_mean(b,:,s,:);
        print_med(1:n_file,1:canali) = ris_med(b,:,s,:);

        for i=1:1:8     %divisione in emisferi
            print_arv_DX(:,i) =  print_arv(:,ordine_DX(i));
            print_arv_SX(:,i) =  print_arv(:,ordine_SX(i));
            print_mean_DX(:,i) =  print_mean(:,ordine_DX(i));
            print_mean_SX(:,i) =  print_mean(:,ordine_SX(i));
            print_med_DX(:,i) =  print_med(:,ordine_DX(i));
            print_med_SX(:,i) =  print_med(:,ordine_SX(i));
        end

        %stampa delle matrici su file .csv (separatore ,)
        %emisfero DX
        csvwrite([PathSave '\SUBJ_x_CH_DX\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_arv_DX);
        csvwrite([PathSave '\SUBJ_x_CH_DX\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_mean_DX);
        csvwrite([PathSave '\SUBJ_x_CH_DX\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_med_DX);
        %emisfero SX
        csvwrite([PathSave '\SUBJ_x_CH_SX\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_arv_SX);
        csvwrite([PathSave '\SUBJ_x_CH_SX\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_mean_SX);
        csvwrite([PathSave '\SUBJ_x_CH_SX\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_sezione_' tasks{s} '.csv'],print_med_SX);
    end       
    end

    for b=1:1:n_band
    for c=1:1:canali

        %cancello le variabili per evitare problemi di dimensione delle matrici dell'operazione precedente
        clear print_arv print_mean print_med print_psd print_max_psd print_X_XII_psd

        %blocco l'incide della banda e quello del canale per salvare soggetti x task 
        print_arv(1:n_file,1:n_sez) = ris_arv(b,:,:,c);
        print_mean(1:n_file,1:n_sez) = ris_mean(b,:,:,c);
        print_med(1:n_file,1:n_sez) = ris_med(b,:,:,c);
        print_psd(1:n_file,1:n_sez) = ris_density(b,:,:,c);
        print_max_psd(1:n_file,1:n_sez) = ris_max_density(b,:,:,c);
        print_X_XII_psd(1:n_file,1:n_sez) = ris_X_XII_density(b,:,:,c);
        
        %stampa delle matrici su file .csv (separatore ,)
        csvwrite([PathSave '\SUBJ_x_TSK\ARV_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_arv);
        csvwrite([PathSave '\SUBJ_x_TSK\MEAN_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_mean);
        csvwrite([PathSave '\SUBJ_x_TSK\MED_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_med);
        csvwrite([PathSave '\SUBJ_x_TSK\PSD_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_psd);
        csvwrite([PathSave '\SUBJ_x_TSK\MAX_PSD_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_max_psd);
        csvwrite([PathSave '\SUBJ_x_TSK\X_XII_PSD_(' num2str(fr_inf(b)) '-' num2str(fr_sup(b)) ')_canale_' num2str(c) '.csv'],print_X_XII_psd);
        
    end       
    end

end


save([PathSave '\info.mat'],'versione','canali','files_original','fr_inf','fr_sup','dim_sec' ,'freq_camp','ordine_DX','ordine_SX','n_band','ris_arv','ris_devstd','ris_mean','ris_med','ris_mean_round','ris_max_density','ris_density', 'ris_X_XII_density');

msgbox('>> >> Operazione Completata << <<','Operazione completata','help');



%Versione 4.2: 

%Versione 4.1.2: ADDON salvataggio file con info sull'analisi effettuata

%Versione 4.1.1: FIX mediana del file soggetti x canali

