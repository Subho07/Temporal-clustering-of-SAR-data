clc; clear all; close all;

[fileName,FolderName] = uigetfile('*.*', 'Path selection Time 1');
cd(FolderName);
C = strsplit(FolderName,'\');
 
config_ID = fopen(strcat(FolderName,'\','config.txt'),'rb');
tline = fgetl(config_ID);
tline = fgetl(config_ID);
b = str2num(tline); %row
tline = fgetl(config_ID);
tline = fgetl(config_ID);
tline = fgetl(config_ID);
a = str2num(tline); %column
 
 
cd(FolderName)
fileList = dir('*.bin');
 
folderName = strcat(FolderName,'\','T11.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T11 = fread(fileID,[a b],'float32');
T11 = T11';

folderName = strcat(FolderName,'\','T12_imag.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T12_imag = fread(fileID,[a b],'float32');
T12_imag = T12_imag';

folderName = strcat(FolderName,'\','T12_real.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T12_real = fread(fileID,[a b],'float32');
T12_real = T12_real';

T12 = complex(T12_real,T12_imag);

folderName = strcat(FolderName,'\','T13_imag.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T13_imag = fread(fileID,[a b],'float32');
T13_imag = T13_imag';

folderName = strcat(FolderName,'\','T13_real.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T13_real = fread(fileID,[a b],'float32');
T13_real = T13_real';

T13 = complex(T13_real,T13_imag);
 
folderName = strcat(FolderName,'\','T22.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T22 = fread(fileID,[a b],'float32');
T22 = T22';

folderName = strcat(FolderName,'\','T23_imag.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T23_imag = fread(fileID,[a b],'float32');
T23_imag = T23_imag';

folderName = strcat(FolderName,'\','T23_real.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T23_real = fread(fileID,[a b],'float32');
T23_real = T23_real';

T23 = complex(T23_real,T23_imag);

folderName = strcat(FolderName,'\','T33.bin');
disp(folderName);
fileID = fopen(folderName,'rb');
T33 = fread(fileID,[a b],'float32');
T33 = T33';
 
theta = zeros(b,a);
ent = zeros(b,a);
alpha = zeros(b,a);

B1 = zeros(b,a);
B2 = zeros(b,a);
B3 = zeros(b,a);
[nrows,ncols]= size(T11);
 %%
 %% for window processing

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window
fprintf('Window size: %d \n',ii);

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column
% Starting row and column fixed by the size of the patch extracted from the image of 21/10/1999

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopj= nrows-inci; % Stop row for window processing
stopi= ncols-incj; % Stop column for window processing

t = cputime;

for ii=startj:stopj
    for jj=starti:stopi

        t11s = mean2(T11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t12s = mean2(T12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t13s = mean2(T13(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        
        t21s = conj(mean2(T12(ii-inci:ii+inci,jj-incj:jj+incj)));%i sample
        t22s = mean2(T22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        t23s = mean2(T23(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        t31s = conj(mean2(T13(ii-inci:ii+inci,jj-incj:jj+incj)));%i sample
        t32s = conj(mean2(T23(ii-inci:ii+inci,jj-incj:jj+incj)));%i sample
        t33s = mean2(T33(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        %Coherency matrix
        T_T1 = [t11s t12s t13s; t21s t22s t23s; t31s t32s t33s];
        
        m1 = real(sqrt(1-(27*(det(T_T1)./(trace(T_T1).^3))))); % DOP Barakat

        span = t11s + t22s + t33s;
        h = (t11s - t22s - t33s);
        g = (t22s + t33s);
        val = (m1.*span.*h)./(t11s.*g+m1.^2.*span.^2);
        the = atand(val);
        theta(ii,jj) = 2*the;
        
        
        %entropy
        [evec_v, eval] = eig(T_T1);
        
       
        eval_diag = (sort(diag(eval)))';
        
        if (eval_diag(1) < 0)
            eval_diag(1) = 0;
        end
        
        if (eval_diag(2) < 0)
            eval_diag(2) = 0;
        end
        
        if (eval_diag(3) < 0)
            eval_diag(3) = 0;
        end
        
        %Lambda 1
        eval_norm1 = (eval_diag(3))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm1(eval_norm1 < 0) = 0;
        eval_norm1(eval_norm1 > 1) = 1;
        
        B1(ii,jj) = eval_norm1;
        
        %Lambda 2
        eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
        B2(ii,jj) = eval_norm2;
        
        %Lambda 3
        eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm3(eval_norm3 < 0) = 0;
        eval_norm3(eval_norm3 > 1) = 1;
        
        B3(ii,jj) = eval_norm3;
        
        %Entropy
        ent(ii,jj) = 1-(-eval_norm1*log10(eval_norm1)./log10(3) - ...
            eval_norm2*log10(eval_norm2)./log10(3) - ...
            eval_norm3*log10(eval_norm3)./log10(3));
    end
    fprintf('Column: %d \n',ii);
end

%%
min_th = min(min(theta));
max_th = max(max(theta));
min_th1 = min(min(ent));
max_th1 = max(max(ent));
%%
 
if ~exist(strcat(FolderName,'\','Theta_SAR_clustering'), 'dir')
       mkdir('Theta_SAR_clustering')
       disp('Diresctory created');
end
fclose('all');
Fold = strcat(FolderName,'\','Theta_SAR_clustering\');


hold off
%%
disp('Saving as bin file')
f_name_100 = char(strcat('theta_',date,'.bin'));
fileandpath_100=f_name_100;
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,theta', 'float32');
fclose('all');

f_name_101 = char(strcat('Ent_',date,'.bin'));
fileandpath_101=f_name_101;
fid_101 = fopen(fileandpath_101,'wb');
fwrite(fid_101,ent', 'float32');
fclose('all');

cd(Fold)
%%
clc;
% Threshold for whole area
disp('Saving theta entropy threshold image')
fig2 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);

thre_mat_TM = ones(b,a);

for i = 1:b
    for j = 1:a
        if ent(i,j)>0.5 && ent(i,j) <= 1
            if theta(i,j) >=-90 && theta(i,j) <-10
                thre_mat_TM(i,j) = 1;
            end
			
            if theta(i,j)>= -10 && theta(i,j)<0
                thre_mat_TM(i,j) = 4;
            end
            if theta(i,j) >=0 && theta(i,j) <20
                thre_mat_TM(i,j) = 7;
            end

            if theta(i,j) >=20 && theta(i,j) <=90
                thre_mat_TM(i,j) = 10;
            end
        end
        if ent(i,j)>0.3 && ent(i,j) <= 0.5
             if theta(i,j) >=-90 && theta(i,j) <-10
                thre_mat_TM(i,j) = 2;
            end

            if theta(i,j) >= -10 && theta(i,j)<0
                thre_mat_TM(i,j) = 5;
            end
            if theta(i,j) >=0 && theta(i,j) <20
                thre_mat_TM(i,j) = 8;
            end

            if theta(i,j) >=20 && theta(i,j) <=90
                thre_mat_TM(i,j) = 11;
            end
        end
        
        if ent(i,j)>=0 && ent(i,j) <= 0.3
             if theta(i,j) >= -90 && theta(i,j) <-10
                thre_mat_TM(i,j) = 3;
            end

            if theta(i,j) >= -10 && theta(i,j)<0
                thre_mat_TM(i,j) = 6;
            end
            if theta(i,j) >= 0 && theta(i,j) <20
                thre_mat_TM(i,j) = 9;
            end

            if theta(i,j) >= 20 && theta(i,j) <=90
                thre_mat_TM(i,j) = 12;
            end
        end
    end
end


set(gca,'FontSize',20)
imagesc(thre_mat_new)
axis('image');
axis off;

mymap = [
    1,0.69,0.69;%1
    0.97, 0.24, 0.24;%2
    0.77, 0.19, 0.19;%3
    0.51,1,0.84;%4
    0.01, 1, 0.67;%5
    0.50, 0.99, 0.51;%6
    0.45, 0.89, 0.46;%7
    0.35, 0.69, 0.58;%8
    0.11,0.44,0.23;%9
    0.51,0.58,1;%10
    0.2,0.32,0.97;%11
    0,0.15,0.96;%12
   ];

colormap(mymap);
numcolors = 12;
n = 1;
colorbar('YTick',[n+0.5*(numcolors-n)/numcolors:(numcolors-n)/numcolors:numcolors],'YTickLabel',{'Z1','Z2','Z3','Z4','Z5','Z6','Z7','Z8','Z9','Z10','Z11','Z12'}, 'YLim', [1 numcolors],'FontSize', 20,'TickLength',0)

file2 =  char(strcat('theta_Ent_threshold_',date));
saveas(fig2,file2,'png');
saveas(fig2,file2,'fig');

f_name_surface = strcat(['threshold_thetaFP_HFP','.bin']);
fileandpath_surface=strcat([FolderName f_name_surface]);
fid_01 = fopen(fileandpath_surface,'wb');
fwrite(fid_01,thre_mat_TM', 'float32');
fclose('all');
