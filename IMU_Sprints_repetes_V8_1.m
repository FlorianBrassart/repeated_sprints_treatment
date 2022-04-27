%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------Programme D'analyse et de génération de rapport du test des sprints répétés-----%
%-------------Développé par Florian Brassart le 12/01/2021----------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Développé pour des fichiers de données en .json
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Modification le 18/05/2021
%Modification des graphs (vitesse par rapport à distance) et ajout d'un
%export image dans le rapport excel directement par la fonction
%"xlsputimage.m"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fonction appelée : 
% read_IMU_JSON.m ou read_IMU_v2
% xlsputimage.m

    % Fichier excel appelé pour générer le rapport
% Model_Rapport_sprint_repetes.xlsx
    % Nécessité de 


% Version Matlab conseillée : 2018a ou plus
% Format de Fichiers de données : .json ou .bin


% Le rangement des Données

% exemple de chemain d'accès DATA_Triées\001_014_CL\Sprints_repetes\IMU

% Dans cet exemple, le dossier à appeler au début du programme est "DATA_Triées"
% Celui ci contient la liste des dossiers només avec l'identifiant du Sujet
% ("001_014_CL") dans l'exemple.

% Chaque sous dossiers contient un dossier comportant le nom du test,
% celui ci doit imperativement avoir l'ortographe suivante
% "Sprints_repetes"

% De plus il contient un fichier text si besoin avec les informations de
% taille de roue en pouce, de frequence d'aquisition et d'angle de
% carrossage (chaque ligne étant les données de chaque joueur avec un
% espace entre chaque valeurs.

% Dans le dossier sprints repétés il y a un autre sous dossier comportant
% le nom que vous voulez (ici : " IMU", dans l'exemple) qui contient les
% fichiers .json obtenus lors de l'aquisition. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear 
close all
clc
rr=dir([cd '/*.mat']);
if ~isempty(rr)
load('result.mat')
end

%% Calcule orientation fauteuil à partir de données IMU à l'avant du fauteuil.
folder_G=uigetdir(path,'choisir le dossier contenant les donnees'); % Dossier contenant toutes les datas triées en sous dossiers
listing_G=dir(folder_G);
dirData = dir([folder_G '/*.txt']);      %# Get the data for the current directory
file_info = {dirData.name}';  %'# Get a list of the files
listing_G=listing_G([listing_G.isdir] & ~strcmp({listing_G.name},'.') & ~strcmp({listing_G.name},'..'));
if~isempty(file_info)
Infos_roue=load(strcat(folder_G,'\',cell2mat(file_info(1,1))));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% seuil=10; % seuil detection sprint 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               
%% ----- Sélection des fichiers
% boulce par sujet 
for g=1:size(listing_G)
    clear nom_dos listing
    listing=dir(strcat(folder_G,'\',listing_G(g,1).name,'\'));% sous dossiers contenant des sous dossier de tests ou de sujet
    path_g=strcat(folder_G,'\',listing_G(g,1).name,'\');
    listing=listing([listing.isdir] & ~strcmp({listing.name},'.') & ~strcmp({listing.name},'..'));
    nom_dos={listing_G(g,1).name};
    disp(nom_dos)
    % Infos
    if~isempty(file_info)
    rayon_roue = Infos_roue(g,1)/78.74;
    FrHZ= Infos_roue(g,2);
    camber_angle = Infos_roue(g,3);
    M= Infos_roue(g,4);
    else
    % Valeur par défaut des caractéristique fauteuil. 
    FrHZ=128; % fréquence d'aquisition des IMUs
    rayon_roue = 0.3304;
    % rayon_roue = 0.3175;
    camber_angle = 17; % angle de carossage des roues du fauteuil en degre
    dist_roue = 0.80;
    M=80;                            
    end
    nom={};
    %% Boucle par dossier test 
    for f = 1: size(listing,1)
        clear nom True3 condition condition1 tf 
%         True3= '100_m_fllorian_29';
        True3= 'Sprints_repetes';
        True4='Sprints_repetes_r';
        nom={listing(f,1).name};
        disp(strcat(nom_dos,'_',nom))
        condition = char(nom_dos);
        condition1 = char(nom);
        tf = strcmp(condition1,True3);
        tf2= strcmp(condition1,True4);
        if tf2==1
            tf=1;
        end
        switch tf
            case 1
                %% Création des chemins d'accés 
                clear pathSW path infos fileList fileListSW cell DATAR1 DATAR2 DATAF DATAF_av DATASW1 DATASW2 start RECAP_sprint
                pathSW=strcat(listing(f,1).folder,'\',listing(f,1).name);
                infos=dir(pathSW);infos=infos([infos.isdir] & ~strcmp({infos.name},'.') & ~strcmp({infos.name},'..'));
                path =strcat(pathSW,'\',infos(end,1).name);
                   % Importation des données de temps de sprints des cellules photoelectriques 
                dirData_tps = dir([pathSW '/*.txt']);      %# Get the data for the current directory
                file_tps = {dirData_tps.name}';  %'# Get a list of the files
                %     file_tps = getAllFiles(pathSW, '*.txt', 0);
                if~isempty(file_tps)
                    Infos_tps=load(strcat(pathSW,'\',cell2mat(file_tps(1,1))));
                end

                %% Ouverture des fichiers IMU 
                dirData_J = dir([path '/*.json']);      %# Get the data for the current directory
                fileList = {dirData_J.name}';  %'# Get a list of the files
                FF=2;
                num_file=5;
                if isempty(fileList)
                    dirData_bin = dir([path '/*.bin']);      %# Get the data for the current directory
                    fileList = {dirData_bin.name}';  %'# Get a list of the files
                    clear FF num_file
                    FF=1;
                    num_file=4;
                end
                if length(fileList)>=num_file-1
%                    

                    switch FF
                        case 2
                            cc1=contains(fileList,'Droite');
                            cc2=contains(fileList,'Gauche');
                            cc3=contains(fileList,'Cale');
                            file1=strcat(path,'\',char(fileList(cc1)));
                            file2=strcat(path,'\',char(fileList(cc2)));
                            file3=strcat(path,'\',char(fileList(cc3)));
%                             
                            if length(fileList)==num_file
                                cc4=contains(fileList,'Torse');
                                file4=strcat(path,'\',char(fileList(cc4)));

                            else
                                file4=strcat(path,'\',cell2mat(fileList(3,1))); % mise en lien du chemin d'accés + le nom du fichier 2
                            end

                            if isempty(file_info)
%                                 file_infos=strcat(path,'\',cell2mat(fileList(end,1))); % mise en lien du chemin d'accés + le nom du fichier 2
                                 cc5=contains(fileList,'information');
                                file_infos=strcat(path,'\',char(fileList(cc5)));
                                INFOS=jsondecode(fileread(file_infos));
                                rayon_roue=(str2num(INFOS.taille_roues)/39.37)/2;
                                camber_angle=str2num(INFOS.angle_carossage);
%                                 dist_roue=str2num(INFOS.largeur_fauteuil_sol);
                            end
                            DATAR1_struct = jsondecode(fileread(file1));
                            DATAR2_struct = jsondecode(fileread(file2));
                            DATAF_av_struct = jsondecode(fileread(file3));
                            DATADos_struct = jsondecode(fileread(file4));

                            DATAR1 = readIMU_JSON(file1,FrHZ,rayon_roue);% DATAR1([1:10 end-10:end],2:end)=0; % fonction permettant de récupérer les données des fichiers .bin sous forme de matrice en argument il faut renseigner (chemain d'accès + nom du fichier, la fréquence à laquelle les données ont étaient enregistrées, la fréquence de coupure pour le filtre lowpass utilisé.
                            DATAR2 = readIMU_JSON(file2,FrHZ,rayon_roue);%DATAR2([1:10 end-10:end],2:end)=0;
                            DATAF_av= readIMU_JSON(file3,FrHZ,rayon_roue); %DATAF_av([1:10 end-10:end],2:end)=0;% IMU placée sur le cadre
                            DATADos= readIMU_JSON(file4,FrHZ,rayon_roue); %DATARDos([1:10 end-10:end],2:end)=0;% IMU placée sur le dos
                        case 1
                            cc1=contains(fileList,'RD');
                            cc2=contains(fileList,'RG');
                            cc3=contains(fileList,'Cadre');
                            file1=strcat(path,'\',char(fileList(cc1)));
                            file2=strcat(path,'\',char(fileList(cc2)));
                            file3=strcat(path,'\',char(fileList(cc3)));
                            if length(fileList)==num_file
                                cc4=contains(fileList,'Tronc');
                                file4=strcat(path,'\',char(fileList(cc4)));
                            else
                                file4=strcat(path,'\',cell2mat(fileList(3,1))); % mise en lien du chemin d'accés + le nom du fichier 2
                            end
                            DATAR1 = readIMU_v2(file1,FrHZ,rayon_roue);
                            DATAR2 = readIMU_v2(file2,FrHZ,rayon_roue);
                            DATAF_av= readIMU_v2(file3,FrHZ,rayon_roue);
                            offset=mean(DATAF_av(10:2*FrHZ,2));
                            DATAF_av(:,2)=DATAF_av(:,2)-offset;
                            DATADos= readIMU_v2(file4,FrHZ,rayon_roue);
                    end
                    if mean(DATAR1(:,11))<0 % détection du sens de rotation des roues pour les remmetres dans le bon sens
                        DATAR1(:,2:end)=-DATAR1(:,2:end);
                    end
                    if mean(DATAR2(:,11))<0
                        DATAR2(:,2:end)=-DATAR2(:,2:end);
                    end
                    % Si les longueurs des matrices des différents capteurs ne sont pas les
                    % même => cropage en fonction de la matrice la plus courte.
                    clear lgt zz
                    lgt(1,1)=length(DATAR1);
                    lgt(2,1)=length(DATAR2);
                    lgt(3,1)=length(DATAF_av);
                    lgt(4,1)=length(DATADos);
                    [zz, ~]=find(lgt==min(lgt));
                    zz=zz(1,1);
                    clear DATA
                    DATA= DATAR1(1:lgt(zz,1),:);
                    clear DATAR1
                    DATAR1= DATA;
                    clear DATA
                    DATA= DATAR2(1:lgt(zz,1),:);
                    clear DATAR2
                    DATAR2 = DATA;
                    clear DATA
                    DATA= DATAF_av(1:lgt(zz,1),:);
                    clear DATAF_av
                    DATAF_av=DATA;
                    clear DATA
                    DATA= DATADos(1:lgt(zz,1),:);
                    clear DATADos
                    DATADos=DATA;
                    clear DATA
                    sizeDATA=size(DATADos);
                    %% Calcul de la vitesse en prenant en compte l'angle de carrossage d'après fuss et al 2012
                    clear  dif signe Mwz Mwxy Mwz1 Mwxy1 tan0 
                    dif=DATAR2(:,7)-DATAR1(:,7);
                    dif(dif==0)=0.0001;                              
                    signe=dif(:,1)./abs(dif(:,1));
                    Mwz= deg2rad(DATAR1(:,7));
                    Mwxy =deg2rad(sqrt((DATAR1(:,5)).^2+(DATAR1(:,6)).^2));
                    tan0=deg2rad(tan(camber_angle));
                    Mwz1= deg2rad(DATAR2(:,7));
                    Mwxy1=deg2rad(sqrt((DATAR2(:,5)).^2+(DATAR2(:,6)).^2));
                    DATAR1(:,13)=(Mwz-(Mwxy.*tan0.*signe)).*rayon_roue.*3.6;
                    DATAR2(:,13)=(Mwz1-(Mwxy1.*tan0.*(-signe))).*rayon_roue.*3.6;                          
                    %% Filtre DATA IMU
                    Wn=(10/FrHZ);
                    [af,bf] = butter(2,Wn,'low');
                    DATAR1(:,2:end) = filtfilt(af,bf,DATAR1(:,2:end));
                    DATAR2(:,2:end) = filtfilt(af,bf,DATAR2(:,2:end));
                    DATADos(:,2:end) = filtfilt(af,bf,DATADos(:,2:end));
                    DATAF_av(:,2:end) = filtfilt(af,bf,DATAF_av(:,2:end));
                    %% Calcul de la vitesse et de l'accélération moyenne
                    clear Force Puissance DATAmean1 Accmean meandist_tot meandist_R1 meandist_R2 Force_eff

                    DATAmean1=(DATAR1(:,13)+DATAR2(:,13))/2;
                    DATAmean2=DATAmean1(1,1)-(DATAmean1(2,1)-DATAmean1(1,1));
                    DATAmean2=[DATAmean2 ;DATAmean1];
                    Accmean = diff(DATAmean2/3.6)*FrHZ; clear DATAmean2
                    meandist_tot=cumtrapz(DATAmean1/3.6)/FrHZ;
                    meandist_R1=cumtrapz(DATAR1(:,13)/3.6)/FrHZ;
                    meandist_R2=cumtrapz(DATAR2(:,13)/3.6)/FrHZ;
                    Force=Accmean*M;
                    Force_eff=Force;
                    Force_eff(Force_eff<0)=0;
                    Puissance=Force_eff.*(DATAmean1/3.6);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %% -_-_-_-_-_-_-_-_-_-_-Sprint Répétés-_-_-_-_-_-_-_-_-_-_-_-_- %%
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    clear IntvitR IntvitR2 IntaccR IntaccR2 Introtvit Inttps synthese_lineaire synthese_temps synthese_Rotation count count1 count2 Picrot_acc Meanrotacc MeanVit Meanacc Pic_vitR Pic_vit MeanVit2 Meanacc2 Pic_vitR2 Pic_vit2 synthese_MAV synthese_MAR Data_vitesse_events_AV Data_vitesse_events
                    % Détection des début de sprint AUTOMATIQUE
                    clear Event Events start start1 meandist DATAmean
                    seuil = max(DATAmean1)/2;
                    DATAF_av(:,2:4)=DATAF_av(:,2:4)*9.81;
                    start1 = find(DATAmean1>seuil);
                    start1 = start1(1,1);
                    while DATAmean1(start1,1)>0.21
                        start1=start1-1;
                    end
                    start = start1;
                    Event=[0 0];
                    Events=0;


                    count=0;
                    count2=0;
                    count1 = 2;
                    for e=start:size(DATAmean1)
                        if DATAmean1(e,1) > seuil
                            count = count+1;
                            Event(count,:) = [DATAmean1(e,1) e];
                        end
                    end
                    clear count
                    count=1;
                    for b=2:size(Event)-1
                        clear Event1 Ind_event
                        if Event(b,2)-Event(b-1,2)>7*FrHZ 
                            count=count+1;
                            Event1 = Event(b,2);
                            while DATAmean1(Event1,1)>0.21
                                Event1=Event1-1;
                            end
                           Events(count,1) = Event1;
                        end
                    end

                    Events(1,1)=start;
                    Events(end+1,1)=length(DATAmean1);
                    meandist(1,1:size(Events)-1)={0};
                    col{1,1}={'[0,0,1]'};
                    col{2,1}={'[0,0.5,0]'};
                    col{3,1}={'[1,0,0]'};
                    col{4,1}={'[0,1,1]'};
                    col{5,1}={'[1,0,1]'};
                    col{6,1}={'[1,0.5,0]'};
                    
%                     FIG_DECOLLAGE = figure('Units','Normalized','Position',[0.1 0.1 0.7 0.7],'Color',[0.76 0.64 0.58],'NumberTitle','off','Name','Vérification du Point de départ du Déplacement (StartDep)','MenuBar','none');movegui('center');
%                         FIGAXIS = axes('Units','Normalized','Color',[1 0.9686 0.9216],'Position',[0.1 0.1 0.65 0.8],'XColor',[.19 .19 .19],'YColor',[.19 .19 .19]);
%                         Curve = plot(DATAmean1(:,1),'o','MarkerFaceColor','k','MarkerSize',2); hold on; ZOOM = 0;grid on
%                         Marker = plot(Events(1:end-1,1),DATAmean1(Events(1:end-1,1)),'o','color',[1 0.3 0],'LineWidth',1.5,'MarkerSize',7);
%                        
%                         BUT_SEUIL = uicontrol('Units','Normalized','BackgroundColor',[0.94 0.67 0.4],'FontSize',12,'Position',[.8 .45 .15 .1],'String','MODIF SEUIL','Style','pushbutton',...
%                           'callback','seuil = inputdlg(''Indiquer le seuil pour la detection auto des pics :'',''Modification du seuil'', [1 50]);seuil = str2num(seuil{:});clear Event Events start start1 meandist DATAmean; DATAF_av(:,2:4)=DATAF_av(:,2:4)*9.81; start1 = find(DATAmean1(:,1)>seuil); start1 = start1(1,1); while DATAmean1(start1,1)>0.21; start1=start1-1;end; start = start1; Event=[0 0]; Events=0; count=0; count2=0; count1 = 2; for e=start:size(DATAmean1); DATAmean(e,1)=(DATAR1(e,13)+DATAR2(e,13))/2; if DATAmean(e,1) > seuil; count = count+1; Event(count,:) = [DATAmean(e,1) e];  end; end; clear count;count=1; for b=2:size(Event)-1; clear Event1 Ind_event;if Event(b,2)-Event(b-1,2)>7*FrHZ;count=count+1; Event1 = Event(b,2);while DATAmean1(Event1,1)>0.21;Event1=Event1-1;end; Events(count,1) = Event1; end; end; Events(1,1)=start; Events(end+1,1)=length(DATAmean1); meandist(1,1:size(Events)-1)={0};delete(Marker); hold on ; Marker = plot(Events(1:end-1,1),DATAmean1(Events(1:end-1,1)),''o'',''color'',[1 0.3 0],''LineWidth'',1.5,''MarkerSize'',7);');
%                        
%                         BUT_VALID = uicontrol('Units','Normalized','BackgroundColor',[0.50 0.87 0.30],'FontSize',12,'Position',[.8 .15 .15 .125],'String','VALIDER','Style','pushbutton',...
%                             'callback','close');
%                         waitfor(FIG_DECOLLAGE);
                        VitesseMoyenne = {}; VitesseDroite = {}; VitesseGauche = {};time={};

                    %% Boucle par evenement
                    clear tpsPR1 tpsPR2 tpsP_mean tpsReR1 tpsReR2 tpsRe_mean tpscycleR1 tpscycleR2 tpscycle_mean Asy cadenceR1 cadenceR2 cadence Acc_push vmoyR1start vmoyR2start vmoystart accR1stab accR1start accR2stab accR2start Acc_push_mean Force_push Acc_push
                    clear END_Events PicvitR vmoyR1stab vmoyR2stab vmoystab vmaxR1start1 vmaxR2start1 Vmoymax1 vmaxR1start2 vmaxR2start2 Vmoymax2 vmaxR1start3 vmaxR2start3 Vmoymax3 vmaxR1stab vmaxR2stab vmaxmoystab Vmax Vmoy accstart Accmax Accmoy accstab loc_max Inttps tpsPR1 tpsPR1 tpsP_mean tpsReR1 tpsReR2 tpsRe_mean tpscycleR1 tpscycleR2 tpscycle_mean Asy cadenceR1 cadenceR2 cadence min_tps
                    for E=1:numel(Events)-1 
                        %% %%%%%%%%%%%%%%%%%%%% VERIFICATION %%%%%%%%%%%%%%%%%%%%%%%%% %%   
                        clear plus Ind_end mat_dist AccEvent StarDep Accmean1
                        StartDep=Events(E,1);
                        % Vérification Graphique du Start des sprints
                        FIG_DECOLLAGE = figure('Units','Normalized','Position',[0.1 0.1 0.7 0.7],'Color',[0.76 0.64 0.58],'NumberTitle','off','Name','Vérification du Point de départ du Déplacement (StartDep)','MenuBar','none');movegui('center');
                        FIGAXIS = axes('Units','Normalized','Color',[1 0.9686 0.9216],'Position',[0.1 0.1 0.65 0.8],'XColor',[.19 .19 .19],'YColor',[.19 .19 .19]);
                        Curve = plot(DATAmean1(:,1),'o','MarkerFaceColor','k','MarkerSize',2); hold on; ZOOM = 0;grid on
                        Marker = plot(Events(E,1),DATAmean1(Events(E,1)),'o','color',[1 0.3 0],'LineWidth',1.5,'MarkerSize',7);
                        title({[cell2mat(nom) '-Essai : ' num2str(E) '/' num2str(length(Events)-1) '}'] ; ' Startsprint ' ; [' Sujet : ' cell2mat(nom_dos)]}) ;
                        BUT_VALUE_Title = uicontrol('Units','Normalized','BackgroundColor',[.60 .60 .60],'FontSize',12,'Position',[.8 .80 .15 .05],'String','Coordonées','Style','text',...
                            'callback','');
                        BUT_VALUE_1 = uicontrol('Units','Normalized','BackgroundColor',[.80 .80 .80],'FontSize',12,'Position',[.8 .75 .15 .05],'String',num2str(3*FrHZ),'Style','text',...
                            'callback','');
                        BUT_VALUE_2 = uicontrol('Units','Normalized','BackgroundColor',[.80 .80 .80],'FontSize',12,'Position',[.8 .70 .15 .05],'String',num2str(DATAmean1(Events(E,1))),'Style','text',...
                            'callback','');
                        BUT_ZOOM = uicontrol('Units','Normalized','BackgroundColor',[.60 .60 .60],'FontSize',12,'Position',[.8 .55 .15 .125],'String','ZOOMER','Style','pushbutton',...
                           'callback','if ZOOM == 0;zoom on;ZOOM = 1;set(BUT_ZOOM,''String'',''DEZOOMER'');elseif ZOOM == 1;zoom off;xlim auto;ylim auto;ZOOM = 0;set(BUT_ZOOM,''String'',''ZOOMER'');end');
                        BUT_MODIF = uicontrol('Units','Normalized','BackgroundColor',[0.94 0.67 0],'FontSize',12,'Position',[.8 .30 .15 .125],'String','MODIFIER','Style','pushbutton',...
                            'callback','[Time,Value]=ginput(1);StartDep = round(Time);delete(Marker);Marker = plot(StartDep,DATAmean1(StartDep),''o'',''color'',[1 0.3 0],''LineWidth'',1.5,''MarkerSize'',7);set(BUT_VALUE_1,''String'',num2str(StartDep));set(BUT_VALUE_2,''String'',num2str(DATAmean1(StartDep)));');
                        
                        %                         if E==1
%                             BUT_SEUIL = uicontrol('Units','Normalized','BackgroundColor',[0.94 0.67 0.4],'FontSize',12,'Position',[.8 .45 .15 .1],'String','MODIF SEUIL','Style','pushbutton',...
%                                 'callback','seuil = inputdlg(''Indiquer le seuil pour la detection auto des pics :'',''Modification du seuil'', [1 50]);seuil = str2num(seuil{:});clear Event Events start start1 meandist DATAmean; DATAF_av(:,2:4)=DATAF_av(:,2:4)*9.81; start1 = find(DATAmean1(:,1)>seuil); start1 = start1(1,1); while DATAmean1(start1,1)>0.21; start1=start1-1;end; start = start1; Event=[0 0]; Events=0; count=0; count2=0; count1 = 2; for e=start:size(DATAmean1); DATAmean(e,1)=(DATAR1(e,13)+DATAR2(e,13))/2; if DATAmean(e,1) > seuil; count = count+1; Event(count,:) = [DATAmean(e,1) e];  end; end; clear count;count=1; for b=2:size(Event)-1; clear Event1 Ind_event;if Event(b,2)-Event(b-1,2)>7*FrHZ;count=count+1; Event1 = Event(b,2);while DATAmean1(Event1,1)>0.21;Event1=Event1-1;end; Events(count,1) = Event1; end; end; Events(1,1)=start; Events(end+1,1)=length(DATAmean1); meandist(1,1:size(Events)-1)={0};');
%                         end
                        BUT_VALID = uicontrol('Units','Normalized','BackgroundColor',[0.50 0.87 0.30],'FontSize',12,'Position',[.8 .15 .15 .125],'String','VALIDER','Style','pushbutton',...
                            'callback','close');
                        waitfor(FIG_DECOLLAGE);

                        Events(E,1)=StartDep;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

                        % boucle par data entre chaque evenements pour mesurer la distance
                        % parcourue depuis le début du sprint 
                        for e=Events(E,1):Events(E+1,1)
                            meandist((e-Events(E,1))+2,E)={(DATAmean1(e,1)*((1/3600)/FrHZ)*1000)+cell2mat(meandist((e-Events(E,1))+1,E))};
                        end
                        mat_dist=cell2mat(meandist(:,E));
                        Ind_end=find(mat_dist>=20);  %%%%%%%% Distance du 20m pour détecter la fin du sprint
                        END_Events(E,1)=Ind_end(1,1)+Events(E,1);

                        %%   % Detection de Pic à partir de l'acceleration et vérification          
                        clear peak_minR1 peak_maxR1 peak_minR2 peak_maxR2 peak_accR1 peak_accR2 DATAfilt_R1 DATAfilt_R2 z x pminr1 pminr2 pmaxr1 pmaxr2 str color VPPP VP z
                        VP(:,1)=DATAR1(:,13);VP(:,2)=DATAR2(:,13);
                        str = '#BCD4E6';
                        color = sscanf(str(2:end),'%2x%2x%2x',[1 3])/255;
                        fig_2 = figure('Units','Normalized','Position',[0.1 0.1 0.7 0.7],'Color',color,'NumberTitle','off','Name','Vérification du Point de départ du Déplacement (StartDep)','MenuBar','none');movegui('center');
                        FIGAXIS = axes('Units','Normalized','Color',[1 0.9686 0.9216],'Position',[0.1 0.1 0.65 0.8],'XColor',[.19 .19 .19],'YColor',[.19 .19 .19]);set(gcf, 'WindowState', 'maximized');
                        % Détection auto des pics
                        z=1;x(1,1)=Events(E,1);x(2,1)=END_Events(E,1); 
                        VPPP=VP(x(1,z):x(2,z),:);
                        %%%%%%%
                        [~,peak_maxR1] = findpeaks(DATAR1(Events(E,1):END_Events(E,1),13), 'MinPeakWidth',12, 'MinPeakDistance', 20, 'MinPeakProminence', 0.5);
                        for i=1:length(peak_maxR1)-1
                            [~,peak_minR1(i+1,1)]=min(VPPP(peak_maxR1(i,1):peak_maxR1(i+1,1),1));
                            peak_minR1(i+1,1)=peak_minR1(i+1,1)+peak_maxR1(i,1)-1;
                        end
                        [~,peak_minR1(1,1)]=min(VPPP(1:peak_maxR1(1),1));
                        [~,peak_minR1(end+1,1)]=min(VPPP(peak_maxR1(end):end,1));
                        peak_minR1(end,1)=peak_minR1(end,1)+peak_maxR1(end)-1;
                        %%%%%%
                        [~,peak_maxR2] = findpeaks(DATAR2(Events(E,1):END_Events(E,1),13), 'MinPeakWidth',12, 'MinPeakDistance', 20, 'MinPeakProminence', 0.5);
                        for i=1:length(peak_maxR2)-1
                            [~,peak_minR2(i+1,1)]=min(VPPP(peak_maxR2(i,1):peak_maxR2(i+1,1),2));
                            peak_minR2(i+1,1)=peak_minR2(i+1,1)+peak_maxR2(i,1)-1;
                        end
                        [~,peak_minR2(1,1)]=min(VPPP(1:peak_maxR2(1),2));
                        [~,peak_minR2(end+1,1)]=min(VPPP(peak_maxR2(end):end,2));
                        peak_minR2(end,1)=peak_minR2(end,1)+peak_maxR2(end)-1;
                        %%%%%%%%%
                        pminr1=peak_minR1;pminr2=peak_minR2;pmaxr1=peak_maxR1;pmaxr2=peak_maxR2;
                        clear peak_minR1 peak_maxR1 peak_minR2 peak_maxR2
                        % Vérification et correction
                        
                        sp1=subplot(2,1,1);
                        plot(VP(x(1,z):x(2,z),1),'r','DisplayName','Roue droite')
                        hold on
                        ppmin1=plot(pminr1(:,1),(VPPP(pminr1(:,1),1)),'o','MarkerFaceColor','m','MarkerSize',10,'DisplayName','right pic min');
                        ppmax1=plot(pmaxr1(:,1),(VPPP(pmaxr1(:,1),1)),'o','MarkerFaceColor','c','MarkerSize',10,'DisplayName','right pic max');
                        legend ('Position',[.83 .60 .06 .06])
                        sp2= subplot(2,1,2);
                        plot(VP(x(1,z):x(2,z),2),'DisplayName','Roue gauche')
                        hold on
                        ppmin2=plot(pminr2(:,1),(VPPP(pminr2(:,1),2)),'o','MarkerFaceColor','r','MarkerSize',10,'DisplayName','left pic min');
                        ppmax2=plot(pmaxr2(:,1),(VPPP(pmaxr2(:,1),2)),'o','MarkerFaceColor','b','MarkerSize',10,'DisplayName','left pic max');   
                        legend ('Position',[.83 .18 .06 .06])
                        BUTADMAXR1 = uicontrol('Units','Normalized','BackgroundColor','c','FontSize',12,'Position',[.0 .835 .12 .08],'String','Ajouter pic max R1','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pmaxr1(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1));pmaxr1(end,1)=pmaxr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pmaxr1(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1));pmaxr1(end,1)=pmaxr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp1;delete(ppmax1);hold on;pmaxr1=sort(pmaxr1) ;sp1=subplot(2,1,1); hold on;ppmax1=plot(pmaxr1(:,1),(VPPP(pmaxr1(:,1),1)),''o'',''MarkerFaceColor'',''c'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .60 .06 .06]);');
                        BUTADMINR1 = uicontrol('Units','Normalized','BackgroundColor','m','FontSize',12,'Position',[.0 .755 .12 .08],'String','Ajouter pic min R1','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pminr1(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1));pminr1(end,1)=pminr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pminr1(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1));pminr1(end,1)=pminr1(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp1;delete(ppmin1);hold on;pminr1=sort(pminr1) ;sp1=subplot(2,1,1); hold on;ppmin1=plot(pminr1(:,1),(VPPP(pminr1(:,1),1)),''o'',''MarkerFaceColor'',''m'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .60 .06 .06]);');
                        BUTSUPMAXR1 = uicontrol('Units','Normalized','BackgroundColor','c','FontSize',12,'Position',[.0 .675 .12 .08],'String','Supprimer pic max R1','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)==max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)==max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pmaxr1(find(pmaxr1==indexpos(i,1)),:)=[];end;sp1;delete(ppmax1);hold on;pmaxr1=sort(pmaxr1) ;sp1=subplot(2,1,1); hold on;ppmax1=plot(pmaxr1(:,1),(VPPP(pmaxr1(:,1),1)),''o'',''MarkerFaceColor'',''c'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .60 .06 .06]);');
                        BUTSUPMINR1 = uicontrol('Units','Normalized','BackgroundColor','m','FontSize',12,'Position',[.0 .595 .12 .08],'String','Supprimer pic min R1','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)==min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),1)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)==min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,1)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pminr1(find(pminr1==indexpos(i,1)),:)=[];end;sp1;delete(ppmin1);hold on;pminr1=sort(pminr1) ;sp1=subplot(2,1,1); hold on;ppmin1=plot(pminr1(:,1),(VPPP(pminr1(:,1),1)),''o'',''MarkerFaceColor'',''m'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .60 .06 .06]);');

                        BUTADMAXR2 = uicontrol('Units','Normalized','BackgroundColor','b','FontSize',12,'Position',[.0 .37 .12 .08],'String','Ajouter pic max R2','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pmaxr2(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2));pmaxr2(end,1)=pmaxr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pmaxr2(end+1,1)]=max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2));pmaxr2(end,1)=pmaxr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp2;delete(ppmax2);hold on;pmaxr2=sort(pmaxr2) ;sp2=subplot(2,1,2); hold on;ppmax2=plot(pmaxr2(:,1),(VPPP(pmaxr2(:,1),2)),''o'',''MarkerFaceColor'',''b'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .60 .06 .06]);');
                        BUTADMINR2 = uicontrol('Units','Normalized','BackgroundColor','r','FontSize',12,'Position',[.0 .29 .12 .08],'String','Ajouter pic min r2','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1;sgtitle(''Clic Gauche sur la courbe pour ajouter des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);[~,pminr2(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2));pminr2(end,1)=pminr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;else;[~,pminr2(end+1,1)]=min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2));pminr2(end,1)=pminr2(end,1)+round(x1(i)-(0.1*FrHZ))-1;end;end;sp2;delete(ppmin2);hold on;pminr2=sort(pminr2) ;sp2=subplot(2,1,2); hold on;ppmin2=plot(pminr2(:,1),(VPPP(pminr2(:,1),2)),''o'',''MarkerFaceColor'',''r'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .60 .06 .06]);');
                        BUTSUPMAXR2 = uicontrol('Units','Normalized','BackgroundColor','b','FontSize',12,'Position',[.0 .21 .12 .08],'String','Supprimer pic max R2','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics max'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)==max(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)==max(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pmaxr2(find(pmaxr2==indexpos(i,1)),:)=[];end;sp2;delete(ppmax2);hold on;pmaxr2=sort(pmaxr2) ;sp2=subplot(2,1,2); hold on;ppmax2=plot(pmaxr2(:,1),(VPPP(pmaxr2(:,1),2)),''o'',''MarkerFaceColor'',''b'',''MarkerSize'',10,''DisplayName'',''right pic max'');legend (''Position'',[.83 .18 .06 .06]);');
                        BUTSUPMINR2 = uicontrol('Units','Normalized','BackgroundColor','r','FontSize',12,'Position',[.0 .13 .12 .08],'String','Supprimer pic min R2','Style','pushbutton',...
                            'callback','MODIF=1;clear x1 y1 indexpos;sgtitle(''Clic Gauche sur la courbe pour supprimer des pics min'');[x1,y1,bouton_souris] = ginput; for i=1:length(x1);if round(x1(i)+(0.1*FrHZ))<length(VPPP);indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)==min(VPPP(round(x1(i)-(0.1*FrHZ)):round(x1(i)+(0.1*FrHZ)),2)));indexpos(i,1)=indexpos(i,1)+round(x1(i)-(0.1*FrHZ))-1;else;indexpos(i,1)=find(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)==min(VPPP(round(x1(i)-(0.1*FrHZ)):end-1,2)));indexpos(i,1)=indexpos+round(x1(i)-(0.1*FrHZ))-1;end;pminr2(find(pminr2==indexpos(i,1)),:)=[];end;sp2;delete(ppmin2);hold on;pminr2=sort(pminr2) ;sp2=subplot(2,1,2); hold on;ppmin2=plot(pminr2(:,1),(VPPP(pminr2(:,1),2)),''o'',''MarkerFaceColor'',''r'',''MarkerSize'',10,''DisplayName'',''right pic min'');legend (''Position'',[.83 .18 .06 .06]);');
                        BUT_VALID = uicontrol('Units','Normalized','BackgroundColor',[0.50 0.87 0.30],'FontSize',12,'Position',[.0 .01 .12 .10],'String','VALIDER','Style','pushbutton',...
                                'callback','close');
                        waitfor(fig_2);

                        peak_minR1(:,1) = pminr1(:,1); 
                        peak_minR1(:,2) = VPPP(peak_minR1(:,1),1);
                        peak_minR2(:,1) = pminr2(:,1); 
                        peak_minR2(:,2) = VPPP(peak_minR2(:,1),2);
                        peak_maxR1(:,1) = pmaxr1(:,1);  
                        peak_maxR1(:,2) = VPPP(peak_maxR1(:,1),1);
                        peak_maxR2(:,1) = pmaxr2(:,1); 
                        peak_maxR2(:,2) = VPPP(peak_maxR2(:,1),2);

                        clear pminr1 pminr2 pmaxr1 pmaxr2  
                        %%    %%%%%%%%%%%%%%%%%%%%%% CALCULS SPRINTS REPETES %%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % sauvegarde vitesses 
                        ccc=0;
                        time_unit = 1/FrHZ;
                        time1 = (time_unit : time_unit: time_unit*length(VPPP))';
                        for cc = x(1,1):x(2,1)
                            ccc = ccc+1;
                        VitesseMoyenne(ccc,E) = {DATAmean1(cc,1)/3.6};
                        VitesseDroite(ccc,E) = {DATAR1(cc,13)/3.6};
                        VitesseGauche(ccc,E) = {DATAR2(cc,13)/3.6};
                        time(ccc,E) = {time1(ccc,1)};
                        end

                        % Calcul de temps 
                        Inttps(E,1) = (END_Events(E,1)-Events(E,1))/FrHZ;
                        if~isempty(file_tps)
                            Inttps_cell(E,1) = Infos_tps(E,1);
                        else
                            Inttps_cell=[]; 
                        end

                        % Calcul vitesse moyenne
                        Vmoy(E,1) = mean(DATAmean1(Events(E,1):END_Events(E,1),1));        
                        % Calcul vitesse pic par interval
                        Pic_vitR(E,1) = max(DATAR1(Events(E,1):END_Events(E,1),13));
                        Pic_vitR(E,2) = (find((DATAR1(Events(E,1):END_Events(E,1),13))== max(DATAR1(Events(E,1):END_Events(E,1),13))))/FrHZ;
                        Pic_vitR(E,3) = max(DATAR2(Events(E,1):END_Events(E,1),13));
                        Pic_vitR(E,4) = (find((DATAR2(Events(E,1):END_Events(E,1),13))== max(DATAR2(Events(E,1):END_Events(E,1),13))))/FrHZ;
                        Vmax(E,1) = max(DATAmean1(Events(E,1):END_Events(E,1),1));
                        % Calcul position de distacne pour atteindre la Vmax
                        clear loc_Ind
                        loc_Ind=find(DATAmean1(Events(E,1):END_Events(E,1),1)== max(DATAmean1(Events(E,1):END_Events(E,1),1)));
                        loc_max(E,1) = mat_dist(loc_Ind);
                        % Calcul Acceleration Moyenne
                        AccEvent=Accmean(Events(E,1):END_Events(E,1),1);
                        if mean(AccEvent(1:3,1))<0
                            AccEvent=-AccEvent;
                        end
                        Accmoy(E,1) = mean(AccEvent);
                        % Calcul Acceleration Max
                        Accmax(E,1) = max(AccEvent);
                        % Calcul entre les 3 premiers pics (phase de demarage) 
                            % vitesse moy 3 1er

                        vmoyR1start(E,1) =mean(DATAR1(Events(E,1):peak_minR1(4,1)+Events(E,1)-1,13));
                        vmoyR2start(E,1) =mean(DATAR2(Events(E,1):peak_minR2(4,1)+Events(E,1)-1,13));
                        vmoystart(E,1) =mean(DATAmean1(Events(E,1):peak_minR2(4,1)+Events(E,1)-1,1));
                            % vitesse max 3 1er
                        vmaxR1start1(E,1)=peak_maxR1(1,1)/FrHZ;
                        vmaxR1start1(E,2)=peak_maxR1(1,2);
                        vmaxR2start1(E,1)=peak_maxR2(1,1)/FrHZ;
                        vmaxR2start1(E,2)=peak_maxR2(1,2);
                        Vmoymax1(E,1)=(peak_maxR1(1,2)+peak_maxR2(1,2))/2;

                        vmaxR1start2(E,1)=peak_maxR1(2,1)/FrHZ;
                        vmaxR1start2(E,2)=peak_maxR1(2,2);
                        vmaxR2start2(E,1)=peak_maxR2(2,1)/FrHZ;
                        vmaxR2start2(E,2)=peak_maxR2(2,2);        
                        Vmoymax2(E,1)=(peak_maxR1(2,2)+peak_maxR2(2,2))/2;

                        vmaxR1start3(E,1)=peak_maxR1(3,1)/FrHZ;
                        vmaxR1start3(E,2)=peak_maxR1(3,2);
                        vmaxR2start3(E,1)=peak_maxR2(3,1)/FrHZ;
                        vmaxR2start3(E,2)=peak_maxR2(3,2);  
                        Vmoymax3(E,1)=(peak_maxR1(3,2)+peak_maxR2(3,2))/2;
                        vmaxstart(E,1) =max(DATAmean1(Events(E,1):peak_minR2(4,1)+Events(E,1)-1,1));
                            % acc moy 3 1er
                        accR1start(E,1) =mean(DATAR1(Events(E,1):peak_minR1(4,1)+Events(E,1)-1,12));
                        accR2start(E,1) =mean(DATAR2(Events(E,1):peak_minR2(4,1)+Events(E,1)-1,12));
                        accstart(E,1) =mean(AccEvent(1:(peak_maxR2(3,1)+peak_maxR1(3,1))/2,1));

                        % Calcul entre les 5 avant-dreniers pics (phase stabilisée)
                            % vitesse moy 5 Der
                        vmoyR1stab(E,1) =mean(DATAR1(peak_minR1(end-5,1)+Events(E,1)-1:peak_minR1(end,1)+Events(E,1)-1,13));
                        vmoyR2stab(E,1) =mean(DATAR2(peak_minR2(end-5,1)+Events(E,1)-1:peak_minR2(end,1)+Events(E,1)-1,13));
                        vmoystab(E,1)=(Vmoymax1(E,1)+Vmoymax2(E,1)+Vmoymax3(E,1))/3;
                            % vitesse max 5 Der
                        vmaxR1stab(E,1) = mean(peak_maxR1(end-5:end,2));
                        vmaxR2stab(E,1) = mean(peak_maxR2(end-5:end,2));
                        vmaxmoystab(E,1) = (vmaxR1stab(E,1)+vmaxR2stab(E,1))/2;
                            % acc moy 5 Der
                        accR1stab(E,1) = mean(DATAR1(peak_minR1(end-5,1)+Events(E,1)-1:peak_minR1(end,1)+Events(E,1)-1,12));
                        accR2stab(E,1) = mean(DATAR2(peak_minR2(end-5,1)+Events(E,1)-1:peak_minR2(end,1)+Events(E,1)-1,12));
                        accstab(E,1)   = mean(AccEvent(peak_minR2(end-5,1):peak_minR2(end,1),1));
                        % Calcul cycles
                        del=length(peak_maxR1)-length(peak_maxR2);
                        if del <= 0
                            sizepeak=length(peak_maxR1);
                        elseif del>0
                            sizepeak=length(peak_maxR2);
                        end
                        Acc_push(1,1)={'Acc Poussée Roue 1'};
                        Acc_push(1,4)={condition};
                        Acc_push(9,1)={'Acc Poussée Roue 2'};
                        Acc_push(17,1)={'temps poussée Roue 1'};
                        Acc_push(17,4)={condition};
                        Acc_push(25,1)={'temps poussée Roue 2'};
                        Acc_push(33,1)={'temps recouvrement Roue 1'};
                        Acc_push(33,4)={condition};
                        Acc_push(41,1)={'temps recouvrement Roue 2'};
                        Acc_push(49,1)={'temps cycle Roue 1'};
                        Acc_push(49,4)={condition};
                        Acc_push(57,1)={'temps cycle Roue 2'};
                        Acc_push(65,1) = {'Force estimée par poussée'};
                        Acc_push(65,4)={condition};
                        Acc_push(73,1) = {'Puissance estimée par poussée'};
                        Acc_push(73,4)={condition};
                        Acc_push(81,1) = {'Vitesse max par poussée roue 1'};
                        Acc_push(81,4)={condition};
                        Acc_push(89,1) = {'Vitesse max poussée roue 2'};
                        Acc_push(97,1) = {'Vitesse min par poussée roue 1 '};
                        Acc_push(97,4)={condition};
                        Acc_push(105,1) = {'Vitesse min par poussée roue 2'};
                        Acc_push(113,4)={condition};
                        Acc_push(113,1) = {'Angle de poussée roue 1'};
                        Acc_push(121,1) = {'Angle de poussée roue 2'};
                        for c=1:sizepeak-1
                            % temps poussée
                            tpsPR1(E,c)={(peak_maxR1(c,1)-peak_minR1(c,1))/FrHZ};
                            tpsPR2(E,c)={(peak_maxR2(c,1)-peak_minR2(c,1))/FrHZ};
                            tpsP_mean(E,c)=(cell2mat(tpsPR1(E,c))+cell2mat(tpsPR2(E,c)))/2;
                            % temps recouvrement
                            tpsReR1(E,c)={(peak_minR1(c+1,1)-peak_maxR1(c,1))/FrHZ};
                            tpsReR2(E,c)={(peak_minR2(c+1,1)-peak_maxR2(c,1))/FrHZ};
                            tpsRe_mean(E,c)=(cell2mat(tpsReR1(E,c))+cell2mat(tpsReR2(E,c)))/2;
                            % temps cycles
                            tpscycleR1(E,c)={(peak_minR1(c+1,1)-peak_minR1(c,1))/FrHZ};
                            tpscycleR2(E,c)={(peak_minR2(c+1,1)-peak_minR2(c,1))/FrHZ};
                            tpscycle_mean(E,c)=(cell2mat(tpscycleR1(E,c))+cell2mat(tpscycleR2(E,c)))/2;
                            % Asymétrie en %age
                                if peak_maxR1(c,2)-peak_maxR2(c,2)>=0
                            Asy(E,c)=(peak_maxR1(c,2)-peak_maxR2(c,2))/peak_maxR1(c,2);
                                else
                            Asy(E,c)=(peak_maxR2(c,2)-peak_maxR1(c,2))/peak_maxR2(c,2);
                                end   
                           % Acceleration par poussée
                           Acc_push(E+1,c)={(((peak_maxR1(c,2)/3.6)-(peak_minR1(c,2)/3.6)))/((peak_maxR1(c,1)-peak_minR1(c,1))/FrHZ)};
                           Acc_push(E+9,c)={(((peak_maxR2(c,2)/3.6)-(peak_minR2(c,2)/3.6)))/((peak_maxR2(c,1)-peak_minR2(c,1))/FrHZ)};
                           Acc_push(E+17,c)=tpsPR1(E,c);
                           Acc_push(E+25,c)=tpsPR2(E,c);
                           Acc_push(E+33,c)=tpsReR1(E,c);
                           Acc_push(E+41,c)=tpsReR2(E,c);
                           Acc_push(E+49,c)=tpscycleR1(E,c);
                           Acc_push(E+57,c)=tpscycleR2(E,c);

                           Acc_push_mean(E,c) = {(cell2mat(Acc_push(E+1,c))+ cell2mat(Acc_push(E+9,c)))/2};
                           Force_push(E,c) = {cell2mat(Acc_push_mean(E,c))*M};
                           Power_push(E,c) = {cell2mat(Force_push (E,c)) * ((peak_maxR1(c,2)/3.6)-(peak_minR1(c,2)/3.6))};
                           
%                            Angle_pushR1(E,c) = {rad2deg((((peak_maxR1(c,2)/3.6)+(peak_minR1(c,2)/3.6))/2)/rayon_roue)*(cell2mat(tpsPR1(E,c)))};
%                            Angle_pushR2(E,c) = {rad2deg((((peak_maxR2(c,2)/3.6)+(peak_minR2(c,2)/3.6))/2)/rayon_roue)*(cell2mat(tpsPR2(E,c)))};

                           Acc_push(E+65,c) = Force_push(E,c);
                           Acc_push(E+73,c) = Power_push(E,c);

                           Acc_push(E+81,c) = {peak_maxR1(c,2)};
                           Acc_push(E+89,c) = {peak_maxR2(c,2)}; 
                           Acc_push(E+97,c) = {peak_minR1(c,2)};
                           Acc_push(E+105,c) = {peak_minR2(c,2)};
%                            Acc_push(E+113,c) = Angle_pushR1(E,c);
%                            Acc_push(E+121,c) = Angle_pushR2(E,c);

                        end    

                        cadenceR1(E,1)=(sizepeak-1)/(peak_minR1(end,1)/FrHZ)*60;
                        cadenceR2(E,1)=(sizepeak-1)/(peak_minR2(end,1)/FrHZ)*60;
                        cadence(E,1)= (cadenceR1(E,1)+cadenceR2(E,1))/2;

                        RECAP_cycle(E,1:3)   = tpsP_mean(E,1:3);
                        RECAP_cycle(E,4)     = mean(tpsP_mean(E,1:3));
                        RECAP_cycle(E,5)     = mean(tpsP_mean(E,end-4:end));
                        RECAP_cycle(E,6)     = mean(tpsP_mean(E,:));
                        RECAP_cycle(E,7:9)   = tpsRe_mean(E,1:3);
                        RECAP_cycle(E,10)    = mean(tpsRe_mean(E,1:3));
                        RECAP_cycle(E,11)    = mean(tpsRe_mean(E,end-4:end));
                        RECAP_cycle(E,12)    = mean(tpsRe_mean(E,:));
                        RECAP_cycle(E,13:15) = tpscycle_mean(E,1:3);
                        RECAP_cycle(E,16)    = mean(tpscycle_mean(E,1:3));
                        RECAP_cycle(E,17)    = mean(tpscycle_mean(E,end-4:end));
                        RECAP_cycle(E,18)    = mean(tpscycle_mean(E,:));
                        RECAP_cycle(E,19:21) = Asy(E,1:3);
                        RECAP_cycle(E,22)    = mean(Asy(E,1:3));
                        RECAP_cycle(E,23)    = mean(Asy(E,end-4:end));
                        RECAP_cycle(E,24)    = cadence(E,1);

                        %%    %%%%%%%%%%%%%%%%%%%%%% FIGURES SPRINTS REPETES %%%%%%%%%%%%%%%%%%%%%%
                        figure(3) % Figure de chaque sprint avec la détection des pics 
                        hold on
                        plot(DATAR1(Events(E,1):END_Events(E,1),13),'k','DisplayName','Right Wheel vit ')% roue droite 
                        plot(peak_minR1(:,1),(peak_minR1(:,2)),'o','MarkerFaceColor','m','MarkerSize',10,'DisplayName','right pic min')
                        plot(peak_maxR1(:,1),(peak_maxR1(:,2)),'o','MarkerFaceColor','c','MarkerSize',10,'DisplayName','right pic max')
                        plot(DATAR2(Events(E,1):END_Events(E,1),13),'b','DisplayName','Left Wheel vit ')% roue gauche
                        plot(peak_minR2(:,1),(peak_minR2(:,2)),'o','MarkerFaceColor','r','MarkerSize',10,'DisplayName','left pic min')
                        plot(peak_maxR2(:,1),(peak_maxR2(:,2)),'o','MarkerFaceColor','b','MarkerSize',10,'DisplayName','left pic max')    
                        legend ('location','northwest')
                        grid on 
                        hold off
                        saveas(figure(3),strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'Sprints_',num2str(E)),'fig')
                        close figure 3
                        clear time_fig count    
                        count=0;
                        % Echelle de Temps en seconde de chaque sprint. 
                        for a=Events(E,1):END_Events(E,1)
                            count=count+1;
                            time_fig(count,1)=count/FrHZ;
                        end    
                        figure(1) % Figure des Sprints Répétés superposées en fonction de la distance
                        set(gcf, 'WindowState', 'maximized');
                        hold on
                        plot(mat_dist(1:END_Events(E,1)-Events(E,1)+1),DATAmean1(Events(E,1):END_Events(E,1),1),'DisplayName',strcat('sprint n°',num2str(E)),'Color',cell2mat(col{E,1}))
%                         plot(mat_dist(1:END_Events(E,1)-Events(E,1)+1),time_fig,'Color',cell2mat(col{E,1}),'HandleVisibility','off')
                        grid on
                    end 
                    legend('location','southeast')
                    set(gca,'fontsize',20);
                    ylabel('vitesse en Km/h','FontSize',25)
                    xlabel('Distance (m)','FontSize',25)
                    xticks(0:1:20)
                    xlim([0 20])
                    hold off 
                    saveas(figure(1),strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'_Superposés'),'fig')
                    saveas(figure(1),strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'_Superposés'),'jpg')
                    imageFileName= strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'_Superposés.jpg');
                    close

                    %% RECAP : 
                    RECAP_sprint(:,1) = vmaxR1start1(:,2);
                    RECAP_sprint(:,2) = vmaxR2start1(:,2);
                    RECAP_sprint(:,3) = Vmoymax1;
                    RECAP_sprint(:,4) = vmaxR1start2(:,2);
                    RECAP_sprint(:,5) = vmaxR2start2(:,2);
                    RECAP_sprint(:,6) = Vmoymax2;
                    RECAP_sprint(:,7) = vmaxR1start3(:,2);
                    RECAP_sprint(:,8) = vmaxR2start3(:,2);
                    RECAP_sprint(:,9) = Vmoymax3;
                    RECAP_sprint(:,10)= vmaxR1stab;
                    RECAP_sprint(:,11)= vmaxR2stab;
                    RECAP_sprint(:,12)= vmaxmoystab;
                    RECAP_sprint(:,13)= vmoystab;
                    RECAP_sprint(:,14)= Vmax;
                    RECAP_sprint(:,15)= Vmoy;
                    RECAP_sprint(:,16)= accstart;
                    RECAP_sprint(:,17)= Accmax;
                    RECAP_sprint(:,18)= Accmoy;
                    RECAP_sprint(:,19)= accstab;
                    RECAP_sprint(:,20)= loc_max;
                    RECAP_sprint(:,22)= vmoystart;
                    RECAP_sprint(:,23:46)=RECAP_cycle;

                     if~isempty(file_tps)
                        RECAP_sprint(:,21)=Inttps_cell;
                        RECAP_sprint(:,47)= Inttps;
                        min_tps=min(Inttps_cell);
                        best_sprint=find(Inttps_cell==min_tps);
                     else
                        RECAP_sprint(:,21)= Inttps;
                        min_tps=min(Inttps);
                        best_sprint=find(Inttps==min_tps);
                     end
                    best_sprint=best_sprint(1,1);
                    recap_best(1,1) = {RECAP_sprint(best_sprint,3)};  % Vmax1
                    recap_best(1,2) = {RECAP_sprint(best_sprint,6)};  % Vmax2
                    recap_best(1,3) = {RECAP_sprint(best_sprint,9)};  % Vmax3
                    recap_best(1,4) = {RECAP_sprint(best_sprint,13)}; % Vmoy (Vmax1,2,3)
                    recap_best(1,5) = {RECAP_sprint(best_sprint,22)}; % Vmoy (Vmax1,2,3)
                    recap_best(1,6) = {RECAP_sprint(best_sprint,14)}; % Vmax
                    recap_best(1,7) = {RECAP_sprint(best_sprint,12)}; % Vmoy Stab
                    recap_best(1,8) = {RECAP_sprint(best_sprint,15)}; % Vmoy
                    recap_best(1,9) = {RECAP_sprint(best_sprint,16)}; % Acc Moy du démarrage
                    recap_best(1,10)= {RECAP_sprint(best_sprint,19)}; % Acc Moy Stab
                    recap_best(1,11)= {RECAP_sprint(best_sprint,18)}; % Acc Moy
                    recap_best(1,12)= {RECAP_sprint(best_sprint,20)}; % Distance/Vmax
                    recap_best(1,13)= {RECAP_sprint(best_sprint,21)}; % Durée
                    recap_best(1,16)= {RECAP_cycle(best_sprint,24)};  % Cadence
                    recap_best(1,17)= {RECAP_cycle(best_sprint,19)};  % Asy1
                    recap_best(1,18)= {RECAP_cycle(best_sprint,20)};  % Asy2
                    recap_best(1,19)= {RECAP_cycle(best_sprint,21)};  % Asy3
                    recap_best(1,20)= {RECAP_cycle(best_sprint,22)};  % Moy (Asy1,2,3)
                    recap_best(1,21)= {RECAP_cycle(best_sprint,23)};  % Asy stab
                    recap_best(1,23)= {RECAP_cycle(best_sprint,4)};   % Moy (Tps P1,P2,P3)
                    recap_best(1,24)= {RECAP_cycle(best_sprint,10)};  % Moy(Tps R1, R2, R3)
                    recap_best(1,25)= {RECAP_cycle(best_sprint,16)};  % Moy (Tps C1, C2, C3)
                    recap_best(1,26)= {RECAP_cycle(best_sprint,5)};   % Tps P Sta 
                    recap_best(1,27)= {RECAP_cycle(best_sprint,11)};  % Tps R sta
                    recap_best(1,28)= {RECAP_cycle(best_sprint,17)};  % Tps Cyc sta
                    recap_best(1,29)= {RECAP_cycle(best_sprint,6)};   % Tps Pou Moy
                    recap_best(1,30)= {RECAP_cycle(best_sprint,12)};  % Tps Rec Moy
                    recap_best(1,31)= {RECAP_cycle(best_sprint,18)};  % Tps Cyc Moy
                    recap_best=recap_best';
                    Best_num(1,1)={['best_sprint : ' num2str(best_sprint)]};
                    nbsprint=length(END_Events);
                    figure(8) % Figure des Sprints Répétés superposées en fonction de la distance
                    set(gcf, 'WindowState', 'maximized');
                    hold on
                    plot(cell2mat(meandist(1:END_Events(best_sprint,1)-Events(best_sprint,1)+1,best_sprint)),DATAmean1(Events(best_sprint,1):END_Events(best_sprint,1),1),'DisplayName',strcat('sprint n°',num2str(best_sprint)),'Color',cell2mat(col{best_sprint,1}))
                    plot(cell2mat(meandist(1:END_Events(nbsprint,1)-Events(nbsprint,1)+1,end)),DATAmean1(Events(nbsprint,1):END_Events(nbsprint,1),1),'DisplayName',strcat('sprint n°',num2str(nbsprint)),'Color',cell2mat(col{nbsprint,1}))
                    legend('location','southeast')
                    set(gca,'fontsize',20);
                    ylabel('vitesse en Km/h','FontSize',25)
                    xlabel('Distance (m)','FontSize',25)
                    xticks(0:1:20)
                    xlim([0 20])
                    grid on 
                    hold off 
                    saveas(figure(8),strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'_best_6'),'fig')
                    saveas(figure(8),strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'_best_6'),'jpg')
                    imageFileName2= strcat(listing(f,1).folder,'\',char(strcat(nom_dos,'_',nom)),'_best_6.jpg');
                    close all
                     %% FIN SWITCH 

                    %% Resultats 
                    OriginDir=cd;
                    % mise des différentes données dans une structure 
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).IMU.Roue_droite = DATAR1;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).IMU.Roue_gauche = DATAR2;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).IMU.cadre = DATAF_av;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).IMU.dos = DATADos;

                    %% sauvgarde résultat evenements Sprints repetes
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).VitesseMoyenne=VitesseMoyenne;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).VitesseDroite=VitesseDroite;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).VitesseGauche=VitesseGauche;


                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).RECAP_sprint=RECAP_sprint;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).RECAP_cycle=RECAP_cycle;
                    result.(strcat('R_',nom_dos{1,1})).(strcat('T_',nom{1,1})).ACC_par_Poussee=Acc_push;

                    clear N FileName_Rapport
                    copyfile([OriginDir '\Model_Rapport_sprints_repetes.xlsx'],listing(f,1).folder);
%                     copyfile([OriginDir '\Model_Rapport_sprint_20M.xlsx'],listing(f,1).folder);        
                    N=split(listing(f,1).folder,'\');
                    FileName_Rapport1 = strcat(N(end,1),'_',listing(f,1).name,'.xlsx');

                    FileName_Rapport=strcat(listing(f,1).folder,'\',cell2mat(FileName_Rapport1));     
                    movefile (strcat(listing(f,1).folder,'\','Model_Rapport_sprints_repetes.xlsx'),FileName_Rapport);   
                    xlswrite(FileName_Rapport,RECAP_sprint,'DATA','I5');
                    xlswrite(FileName_Rapport,recap_best(:,1),'DATA','A4');
                    xlswrite(FileName_Rapport,Best_num(1,1),'DATA','G2');
                    xlswrite(FileName_Rapport,Acc_push,'Recap_poussées');

%                     FileName_Rapport_2=strcat(listing(f,1).folder,'\',cell2mat(N(end,1)),'_',condition1(1:end-7),'de_20m.xlsx');  
%                     movefile (strcat(listing(f,1).folder,'\','Model_Rapport_sprint_20M.xlsx'),FileName_Rapport_2);
%                     xlswrite(FileName_Rapport_2,RECAP_sprint(best_sprint,:),'DATA','N4');
%                     xlswrite(FileName_Rapport_2,Best_num(1,1),'DATA','G2');
                    % Export excel du tableau recapitulatif des poussées uniquement.
                    % Chaque sujet ajoute une fenetre dans le fichier excel. 
                    % dir_acc=cd;
                    % xlswrite(strcat(dir_acc,'\Recap_poussees.xlsx'),Acc_push,condition)
                    imageSize=[482, 326];
                    xlsputimage(FileName_Rapport, imageFileName, 'Rapport Athlète', 'A76', imageSize)
                    xlsputimage(FileName_Rapport, imageFileName2, 'Rapport Athlète', 'A51', imageSize)

                    disp(strcat('fin exportation pour ',nom_dos,'_',nom))
                else
                    disp( ' il est possible qu il n y ai pas assez de capteurs dans le dossier')
                end
            case 0 
                disp( ['Le nom du dossier test est ' condition1 ' et il est différent du nom attendu : ' True3 ' ou ' True4] )
        end
    end
end
% dir_FV='D:\Utilisateurs\Documents\These_Florian\DATA_PARAPERF\IMU';
% xlswrite(strcat(dir_FV,'\Result_FV.xlsx'),DATAFV)
save(strcat(cd,'\','result.mat'),'result');
disp('Tratement terminé')