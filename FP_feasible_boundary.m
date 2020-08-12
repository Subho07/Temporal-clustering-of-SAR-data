% Full / quad pol
gam1 = zeros(1001,1);
ent1 = zeros(1001,1);
B11 = zeros(1001,1);
B21 = zeros(1001,1);
B31 = zeros(1001,1);
B41 = zeros(1001,1);
B51 = zeros(1001,1);
B61 = zeros(1001,1);
B71 = zeros(1001,1);
B81 = zeros(1001,1);
B91 = zeros(1001,1);
k = 1;
for m = 0:0.0005:1
   T = [1,0,0;0,m,0;0,0,m];
   m1 = real(sqrt(1-(27*(det(T)./(trace(T).^3))))); % DOP Barakat
   
    t11s = T(1,1);
    t22s = T(2,2);
    t33s = T(3,3);
   
   span = t11s + t22s + t33s;
   h = (t11s - t22s - t33s);
   g = (t22s + t33s);
        
   val = (m1.*span.*h)./(t11s.*g+m1.^2.*span.^2);
   
   val = atand(val);
   
   gam1(k,1) = 2*val;
   %entropy
        [evec_v, eval] = eig(T);
        
       
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
        
        
        %Lambda 2
        eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
       
        %Lambda 3
        eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm3(eval_norm3 < 0) = 0;
        eval_norm3(eval_norm3 > 1) = 1;
        
        %Alpha 1
        eig_vec_r1 = real(evec_v(1,3));
        eig_vec_c1 = imag(evec_v(1,3));
        
        alpha1 = acos(sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1));
        B41(k,1) = alpha1*180./pi;
        
        %Alpha 2
        eig_vec_r2 = real(evec_v(1,2));
        eig_vec_c2 = imag(evec_v(1,2));
        
        alpha2 = acos(sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2));
        B51(k,1) = alpha2*180./pi;
        
        %Alpha 3
        eig_vec_r3 = real(evec_v(1,1));
        eig_vec_c3 = imag(evec_v(1,1));
        
        alpha3 = acos(sqrt(eig_vec_r3*eig_vec_r3 + eig_vec_c3*eig_vec_c3));
        B61(k,1) = alpha3*180./pi;
        
        %Cloude Alpha
        B71(k,1) = (eval_norm1*alpha1*180./pi + eval_norm2*alpha2*180./pi + ...
            eval_norm3*alpha3*180./pi);
        
        %Entropy
        ent1(k,1) = -eval_norm1*log10(eval_norm1)./log10(3) - ...
            eval_norm2*log10(eval_norm2)./log10(3) - ...
            eval_norm3*log10(eval_norm3)./log10(3);
%         ent1(k,1) = (1-m1);
        data = 1 - m1;
        if data > 0.6058 && data < 0.6062
            disp(t11s);
            disp(t22s);
            disp(t33s);
        end
        k = k + 1;
end

gam1 = gam1(:);
gam1 = deg2rad(gam1);
ent1 = 1 - ent1;
figure('Renderer', 'painters', 'Position', [0 -300 900 900])
polarplot(gam1, ent1, 'LineWidth',2,'Color',[0,0,0])

ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
ax.RLim = [0,1];
ax.ThetaLim = [-90 90];
ax.ThetaTick = [-90 -60 -30 -10 0 20 30 60 90];
ax.FontSize = 20;
ax.LineWidth = 2.5;
ax.Color = 'white';
ax.GridColor = 'black';

hold on

%%
gam2 = zeros(501,1);
ent2 = zeros(501,1);
B12 = zeros(501,1);
B22 = zeros(501,1);
B32 = zeros(501,1);
B42 = zeros(501,1);
B52 = zeros(501,1);
B62 = zeros(501,1);
B72 = zeros(501,1);
B82 = zeros(501,1);
B92 = zeros(501,1);
k = 1;
for m = 0.5:0.001:1
   T = [2*m-1,0,0;0,1,0;0,0,1];
    m1 = real(sqrt(1-(27*(det(T)./(trace(T).^3))))); % DOP Barakat
   
    t11s = T(1,1);
    t22s = T(2,2);
    t33s = T(3,3);
   
   span = t11s + t22s + t33s;
   h = (t11s - t22s - t33s);
   g = (t22s + t33s);
        
   val = (m1.*span.*h)./(t11s.*g+m1.^2.*span.^2);
   
   val = atand(val);
   
   val1 = 2*val;
   if val1 < -90
       val1 = -90;
   end
    gam2(k,1) = val1;

   %entropy
        [evec_v, eval] = eig(T);
        
       
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
        
        
        %Lambda 2
        eval_norm2 = (eval_diag(2))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm2(eval_norm2 < 0) = 0;
        eval_norm2(eval_norm2 > 1) = 1;
        
       
        %Lambda 3
        eval_norm3 = (eval_diag(1))./(eval_diag(1) + eval_diag(2) + eval_diag(3));
        
        eval_norm3(eval_norm3 < 0) = 0;
        eval_norm3(eval_norm3 > 1) = 1;
        
        %Alpha 1
        eig_vec_r1 = real(evec_v(1,3));
        eig_vec_c1 = imag(evec_v(1,3));
        
        alpha1 = acos(sqrt(eig_vec_r1*eig_vec_r1 + eig_vec_c1*eig_vec_c1));
        B42(k,1) = alpha1*180./pi;
        
        %Alpha 2
        eig_vec_r2 = real(evec_v(1,2));
        eig_vec_c2 = imag(evec_v(1,2));
        
        alpha2 = acos(sqrt(eig_vec_r2*eig_vec_r2 + eig_vec_c2*eig_vec_c2));
        B52(k,1) = alpha2*180./pi;
        
        %Alpha 3
        eig_vec_r3 = real(evec_v(1,1));
        eig_vec_c3 = imag(evec_v(1,1));
        
        alpha3 = acos(sqrt(eig_vec_r3*eig_vec_r3 + eig_vec_c3*eig_vec_c3));
        B62(k,1) = alpha3*180./pi;
        
        %Cloude Alpha
        B72(k,1) = (eval_norm1*alpha1*180./pi + eval_norm2*alpha2*180./pi + ...
            eval_norm3*alpha3*180./pi);
        
        %Entropy
        ent2(k,1) = -eval_norm1*log10(eval_norm1)./log10(3) - ...
            eval_norm2*log10(eval_norm2)./log10(3) - ...
            eval_norm3*log10(eval_norm3)./log10(3);
%         ent2(k,1) = (1-m1);
        k = k + 1;
end

gam2 = gam2(:);
gam2 = deg2rad(gam2);
ent2 = 1 - ent2;
polarplot(gam2,ent2,'LineWidth',2,'Color',[0,0,0])


const_gam = deg2rad(-90);
gam5 = repelem(const_gam,100)';
ent5 = linspace(max(ent2),1)';

const_gam1 = deg2rad(20);
gam6 = repelem(const_gam1,100)';
ent6 = linspace(0,1)';

const_gam2 = deg2rad(-10);
gam7 = repelem(const_gam2,100)';
ent7 = linspace(0,1)';

gam8 = repelem(0,100)';
ent8 = linspace(0,1)';

polarplot(gam5,ent5,'LineWidth',2,'Color',[0,0,0])

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
