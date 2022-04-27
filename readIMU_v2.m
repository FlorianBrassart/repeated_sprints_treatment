%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ---------- FONCTION DE LECTURE DES FICHIERS .bin ISSUS ------------------
% ------------------------ DES IMU ATOUTNOVATION --------------------------
% -------------------------------------------------------------------------
% ---------- Développée par Tong Li et Didier PRADON ----------------------
% ---------- le 29/09/2018 ------------------------------------------------
% -------------------------------------------------------------------------
% ---------- Modifiée par Florian -----------------------------------------
% ---------- le 26/02/2019 ------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
% filename = chemin/nom du fichier.bin
% frhz = fréquence d'acquisition
% rayon_roue = pour le calcul de la vitesse en km/h indiquer
% le rayon de la roue en m
% OUPUT
% DATA
% Colonne 1 : temps en sec
% Colonne 2 : accélération en X
% Colonne 3 : accélération en Y
% Colonne 4 : accélération en Z
% Colonne 5 : gyromètre en X
% Colonne 6 : gyromètre en Y
% Colonne 7 : gyromètre en Z
% Colonne 8 : magnétomètre en X
% Colonne 9 : magnétomètre en Y
% Colonne 10 : magnétomètre en Z
% Colonne 11 : vitesse en km/h
% Colonne 12 : accélération linéaire

% Exemple
% filename = 'C:\Users\Kirespi\Documents\MATLAB\fonction ATN capteurs\S2_E6D1056849DF.bin'
% frhz = 128
% rayon_roue = 0.3
% DATA = readIMU_v2(filename,frhz,rayon_roue)
% -------------------------------------------------------------------------
% ---------- Modifiée par Florian le 26/02/2019

% ajout du calcul de l'accélération linéaire à partir de la vitesse (DATA colonne 12)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DATA = readIMU_v2(filename,frhz,rayon_roue)

clear fid data_bin data_bin_brute

%--- Import et lecture des fichier .bin
fid = fopen(filename);
data_bin_brute = fread(fid,[16,inf],'int16','b');
fclose(fid);
data_bin = transpose(data_bin_brute);

%--- Time
clear time_unit time
% [ligne, colone] = size(data_bin);
time_unit = 1/frhz;
time = [time_unit : time_unit: time_unit*(size(data_bin,1))];
time = time';

%--- Stockage des datas dans Matrice
% tps
DATA(:,1) = time;
% accélero
DATA(:,2:4) = data_bin(:,1:3)*(32/2^16);
% gyro
DATA(:,5:7) = data_bin(:,5:7)*(4000/2^16);
% magneto
DATA(:,8:10) = data_bin(:,7:9);

if ~isempty(rayon_roue)
    % calcul de la vitesse de déplacement du fauteuil
    clear vitesse
    vitesse = rayon_roue*(deg2rad(DATA(:,7)))*3.6;
    DATA(:,11) = vitesse;
    
    % calcul accélération linéaire à partir de la vitesse
    clear Acc_frm vit_frm
    vit_frm(1,1) = DATA(1,11)-(DATA(2,11)-DATA(1,11));
    vit_frm(2:size(DATA)+1,1) = DATA(1:end,11);
    vit_frm(end+1,1) = DATA(end,11)+(DATA(end,11)-DATA(end-1,11));
    for i = 2:size(vit_frm,1)-1
        Acc_frm(i-1,1) = ((vit_frm(i+1,1)/3.6)-(vit_frm(i-1,1)/3.6))/(2/frhz);
    end
    DATA(:,12) = Acc_frm;  
end

end

