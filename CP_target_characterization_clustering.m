clc;clear all;close all

[filename, path] = uigetfile('*.*', 'Path selection Time 1');
cd(path);
disp(path);
f0 = fopen([path 'config.txt']);
tmp = fgets(f0);
nrows = sscanf(fgets(f0),'%d');
tmp = fgets(f0);
tmp = fgets(f0);
ncols = sscanf(fgets(f0),'%d');


ep = 0;

%% C2 matrix
f1 = fopen([path 'C11.bin'],'rb');
f2 = fopen([path 'C12_real.bin'],'rb');
f3 = fopen([path 'C12_imag.bin'],'rb');
f4 = fopen([path 'C22.bin'],'rb');

C11 = fread(f1,[ncols nrows],'float32') + ep;
C12 = complex( fread(f2,[ncols nrows],'float32') , fread(f3,[ncols nrows],'float32')) + ep;
C21 = conj(C12);
C22 = fread(f4,[ncols nrows],'float32') + ep;

fclose('all');

psi_in = input('Theta (Orientation): '); % orientation (same as theta in C3-C2_GTLR)

chi_in = input('Tau (Ellipticity: +: LC & -: RC): '); % RC = -theta and LC = +theta

%% for window processing

wsi=input('Window Size: ');
wsj = wsi; % Number of columns in the window

inci=fix(wsi/2); % Up & down movement margin from the central row
incj=fix(wsj/2); % Left & right movement from the central column

starti=fix(wsi/2)+1; % Starting row for window processing
startj=fix(wsj/2)+1; % Starting column for window processing

stopi= nrows-inci; % Stop row for window processing
stopj= ncols-incj; % Stop column for window processing

%% Stokes Parameter

Temp_OC_save = zeros(ncols,nrows);
Temp_SC_save = zeros(ncols,nrows);
Temp_span_save = zeros(ncols,nrows);

theta = zeros(ncols,nrows);
ent = zeros(ncols,nrows);
chi_ang = zeros(ncols,nrows);
alpha_cloude = zeros(ncols,nrows);

B1 = zeros(ncols,nrows);
B2 = zeros(ncols,nrows);
B3 = zeros(ncols,nrows);

t_start = tic;

for ii=startj:stopj
    for jj=starti:stopi
        
        Bin_C11 = mean2(C11(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        Bin_C12 = mean2(C12(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        Bin_C21 = mean2(C21(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        Bin_C22 = mean2(C22(ii-inci:ii+inci,jj-incj:jj+incj));%i sample
        
        c2 = [Bin_C11,Bin_C12;Bin_C21,Bin_C22];
        % Stokes Parameter
        s0 = Bin_C11 + Bin_C22;
        s1 = Bin_C11 - Bin_C22;
        s2 = (Bin_C12 + Bin_C21);
        
        if (chi_in >= 0)
            s3 = (1i.*(Bin_C12 - Bin_C21)); %The sign is according to RC or LC sign !!
        end
        
        if (chi_in < 0)
            s3 = -(1i.*(Bin_C12 - Bin_C21)); %The sign is according to RC or LC sign !!
        end
                
        dop = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0);
        
        SC = ((s0)-(s3))./2;
        Temp_SC_save(ii,jj) = SC;
        OC = ((s0)+(s3))./2;
        Temp_OC_save(ii,jj) = OC;
        span = SC + OC;

        
        h = (OC-SC);
        val = (dop.*span.*h)./(SC*OC+dop.^2.*span.^2);
        the = atand(val);
        theta(ii,jj) = 2*the;
        
        a_cloude = 0.5*atan2d(sqrt(s1^2+s2^2),s3);
        alpha_cloude(ii,jj) = real(a_cloude);
        %entropy
        [evec_v, eval] = eig(c2);
        
       
        eval_diag = (sort(diag(eval)))';
        
        if (eval_diag(1) < 0)
            eval_diag(1) = 0;
        end
        
        if (eval_diag(2) < 0)
            eval_diag(2) = 0;
        end

        %Lambda 1
        eval_norm1 = (eval_diag(2))./(eval_diag(1) + eval_diag(2));
        
        eval_norm1(eval_norm1 < 0) = 0;
        eval_norm1(eval_norm1 > 1) = 1;
        
        B1(ii,jj) = eval_norm1;
        
        %Lambda 2
        eval_norm2 = (eval_diag(1))./(eval_diag(1) + eval_diag(2));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
        B2(ii,jj) = eval_norm2;
        
        %Entropy
        ent(ii,jj) = 1-(-eval_norm1*log10(eval_norm1)./log10(2) - ...
            eval_norm2*log10(eval_norm2)./log10(2));
        
        dop = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0);
        DOCP = (-s3)./(dop.*s0); %changed
        Chi = 0.5.*((180/pi).*asin(DOCP));
        chi_ang(ii,jj) = Chi;
    end
    fprintf('Column: %d \n',ii);
end
%%
min_th = min(min(theta));
max_th = max(max(theta));
min_th1 = min(min(ent));
max_th1 = max(max(ent));
%%
cd(path)
if ~exist(strcat(path,'\','Theta_SAR_HYB_clustering'), 'dir')
       mkdir('Theta_SAR_HYB_clustering')
end
fclose('all');
FolderName = path;
Fold = strcat(FolderName,'\','Theta_SAR_HYB_clustering\');
%%
clc;
disp('Saving theta image')

cd(Fold);
theta_new = theta';

ent_new = ent';

fig1 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
set(gca,'FontSize',20)

imagesc(theta_new);
axis('image');
axis off;

c = jet;
c = flipud(c);
colormap(c);
colorbar('YTick', -90:15:90,'FontSize', 20);
caxis([-90 90]);

C = strsplit(path,'\');

date = 'data';
file1 =  char(strcat('theta_HYB_',date));
saveas(fig1,file1,'png')
saveas(fig1,file1,'fig')
%%
cd(path);
cd(Fold);
disp('Saving as bin file')
f_name_100 = char(strcat('theta_HYB_',date,'.bin'));
fileandpath_100=f_name_100;
fid_100 = fopen(fileandpath_100,'wb');
fwrite(fid_100,theta_new', 'float32');
fclose('all');

f_name_101 = char(strcat('ent_HYB_',date,'.bin'));
fileandpath_101=f_name_101;
fid_101 = fopen(fileandpath_101,'wb');
fwrite(fid_101,ent_new', 'float32');
fclose('all');

cd(Fold)
%%
%H-theta threshold
disp('Saving theta entropy threshold image')
thre_mat = ones(ncols,nrows);


for i = 1:ncols
    for j = 1:nrows
        if ent(i,j)>=0.5 && ent(i,j) <= 1
            if theta(i,j) >-90 && theta(i,j) <-10
                thre_mat(i,j) = 1;
            end

            if theta(i,j)>= -10 && theta(i,j)<0
                thre_mat(i,j) = 4;
            end
            
            if theta(i,j) >=0 && theta(i,j) <20
                thre_mat(i,j) = 7;
            end

            if theta(i,j) >20 && theta(i,j) <90
                thre_mat(i,j) = 10;
            end
        end
        if ent(i,j)>0.3 && ent(i,j) <= 0.5
             if theta(i,j) > -90 && theta(i,j) <-10
                thre_mat(i,j) = 2;
            end

            if theta(i,j) >= -10 && theta(i,j)<0
                thre_mat(i,j) = 5;
            end
            if theta(i,j) >=0 && theta(i,j) <20
                thre_mat(i,j) = 8;
            end

            if theta(i,j) >=20 && theta(i,j) <90
                thre_mat(i,j) = 11;
            end
        end
        
        if ent(i,j)>0 && ent(i,j) <= 0.3
             if theta(i,j) >= -90 && theta(i,j) <-10
                thre_mat(i,j) = 3;
            end

            if theta(i,j) >= -10 && theta(i,j)<0
                thre_mat(i,j) = 6;
            end
            if theta(i,j) >= 0 && theta(i,j) <20
                thre_mat(i,j) = 9;
            end

            if theta(i,j) >= 20 && theta(i,j) <90
                thre_mat(i,j) = 12;
            end
        end
    end
end


fig2 = figure('Renderer', 'painters', 'Position', [0 -300 900 900]);
set(gca,'FontSize',30)
imagesc(thre_mat)
axis('image');
axis off;


mymap = [
    1,1,1;%0
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

file2 =  char(strcat('theta_Ent_threshold_HYB_',date));
saveas(fig2,file2,'png')
saveas(fig2,file2,'fig')