
gam1 = zeros(501,1);
ent1 = zeros(501,1);
chi_in = input('Tau (Ellipticity: +: LC & -: RC): ');
k = 1;
for m = 0:0.001:1/2
        Bin_C11 = (2*m+1)/4;
        Bin_C12 = 1i*((2*m-1)/4);
        Bin_C21 = conj(Bin_C12);
        Bin_C22 = Bin_C11;
        
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
                
        SC = ((s0)-(s3))./2;
        OC = ((s0)+(s3))./2;
        span = SC + OC;
        
        Temp_dop = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0);
        
        h = (SC-OC);
        val = ((Temp_dop*span.*h))./((SC*OC + (Temp_dop.^2)*span.^2));
        
        gamma = atand(val);
        gam1(k,1) = 2*gamma;

        
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
        
        %Lambda 2
        eval_norm2 = (eval_diag(1))./(eval_diag(1) + eval_diag(2));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
        %Entropy
        ent1(k,1) = -eval_norm1*log10(eval_norm1)./log10(2) - ...
            eval_norm2*log10(eval_norm2)./log10(2);
        
        k = k + 1;
end
gam1 = deg2rad(gam1);
ent1 = 1 - ent1;
figure('Renderer', 'painters', 'Position', [0 -300 900 900])
polarplot(gam1,ent1,'LineWidth',2,'Color',[0,0,0])
hold on
%%

gam2 = zeros(501,1);
ent2 = zeros(501,1);
k = 1;
for m = 0:0.001:1/2
        Bin_C11 = (2*m+1)/4;
        Bin_C12 = -1i*((2*m-1)/4);
        Bin_C21 = conj(Bin_C12);
        Bin_C22 = Bin_C11;
        
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
                
        SC = ((s0)-(s3))./2;
        OC = ((s0)+(s3))./2;
        span = SC + OC;
        
        Temp_dop = sqrt((s1).^2 + (s2).^2 + (s3).^2)./(s0);
        
        h = (SC-OC);
        val = ((Temp_dop*span.*h))./((SC*OC + (Temp_dop.^2)*span.^2));
        
        gamma = atand(val);
        gam2(k,1) = 2*gamma;

        
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
        
        %Lambda 2
        eval_norm2 = (eval_diag(1))./(eval_diag(1) + eval_diag(2));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
        %Entropy
        ent2(k,1) = -eval_norm1*log10(eval_norm1)./log10(2) - ...
            eval_norm2*log10(eval_norm2)./log10(2);
        
        k = k + 1;
end

gam2 = deg2rad(gam2);
ent2 = 1 - ent2;
polarplot(gam2,ent2,'LineWidth',2,'Color',[0,0,0])

const_gam1 = deg2rad(20);
gam6 = repelem(const_gam1,100)';
ent6 = linspace(0,1)';

const_gam2 = deg2rad(-10);
gam7 = repelem(const_gam2,100)';
ent7 = linspace(0,1)';

gam8 = repelem(0,100)';
ent8 = linspace(0,1)';

polarplot(gam6,ent6,'k','LineStyle','--','LineWidth',2);
polarplot(gam7,ent7,'k','LineStyle','--','LineWidth',2);
polarplot(gam8,ent8,'k','LineStyle','--','LineWidth',2);

window = 400;
interval_theta = 180/window;
theta_val1 = [-90:interval_theta:90];
theta_val1 = deg2rad(theta_val1);
thre_val1 = zeros((window+1), 1);
thre_val1(:) = 0.5;
thre_val2 = zeros((window+1), 1);
thre_val2(:) = 0.3;


polarplot(theta_val1,thre_val1,'k','LineStyle','--','LineWidth',2);
polarplot(theta_val1,thre_val2,'k','LineStyle','--','LineWidth',2);

ax = gca;
ax.ThetaZeroLocation = 'top';
% d = ax.ThetaDir;
ax.ThetaDir = 'clockwise';
ax.RLim = [0,1];
ax.ThetaLim = [-90 90];
ax.ThetaTick = [-90 -60 -30 -10 0 20 30 60 90];
ax.FontSize = 20;
ax.LineWidth = 2.5;
ax.Color = 'white';
ax.GridColor = 'black';

hold off