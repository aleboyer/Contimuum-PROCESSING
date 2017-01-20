% from frequency spectrum to estimate of SSH horizontal wavenumber spectrum 
% 
% F_eta (k_h) = 2 E_0 b^2 N_0 f / ( ? g^2) ...
%               ?_j [ N_j C_j^5 k_h^2 H(j) / (c_j^2 K_h^2 + f^2)^2  ]
% 
%
addpath /Users/aleboyer/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/MatlabLocal/GarrettMunk3
addpath /Users/aleboyer/ARNAUD/SCRIPPS/IWISE/Matlablocal/matlab_archive/ToolBox-MMP-version1/seawater_ver3_2/
addpath /Users/aleboyer/ARNAUD/TOOLBOXES/Farrar_continuum/GarrettMunk3_ssh/


load('../alb_mat/MD1.mat')
load('../alb_mat/MD1_spec_cline.mat','KEmoor_cline')    
load('../alb_mat/MD2_Continuum_cline_ccw.mat')
load('../alb_mat/MD1_peak_cw_ccw_cline.mat')
load('../alb_mat/GMssh.mat')


%params= Gk91Params;
params= Gm76Params;
kx = logspace(-6,-1,200); 
%kz = logspace(-4,1,200);

ind_moor=1:length(MD);
ind_moor=ind_moor(ind_moor~=619 & ind_moor~=606 & ind_moor~=675);


for i=ind_moor(220:end);
    if (~isempty(peak(i).Continuum)) 
        fprintf('i=%i\n',i)
        for ii=1:length(peak(i).Continuum)
%        for ii=1
            if ~isempty(peak(i).Continuum{ii})
                fprintf('ii=%i\n',ii)
                for iii=1:length(peak(i).Continuum{ii})
%                for iii=1
                    if ~isempty(peak(i).Continuum{ii}{iii})
                    fprintf('iii=%i\n',iii)


% I   = 1  % mooring ID;
% II  = 1  % depth ID;
% III = 1  % time block ID;

N     = MD(i).Nnot*MD(i).Nscale(ii);
om    = logspace(log10(1e-8),log10(N),200);

f     = sw_f(MD(i).lat);
E0    = GMContinuum(i).e0{ii}{iii};
KE    = KEmoor_cline(i).spec{ii}{iii};
freq  = KEmoor_cline(i).freq{ii};
H     = MD(i).waterdepth(ii);
g     = sw_g(MD(i).lat,0);
freq1 = 2*pi*freq./86400;
A     = (N*H/2/2/pi/g)^2;

n=1; % mode number
cg=H/n/pi * (freq1.^2-f^2).^(1/2).*(N^2-freq1.^2).^(3/2)./freq1./(N^2-f^2);
kH=n*pi/H* ((freq1.^2-f^2)./(N^2-freq1.^2)).^(1/2);

K=length(kx);
O=length(om);
F=length(freq);
T=length(E0);

Sx   =zeros(K,T);
So   =zeros(O,T);
SOobs =zeros(F,T);
SKobs =zeros(F,T);
for t=1:T;
    params.E0=E0(t);
    ke=KE(:,t);
    
    Sx(:,t) = GmKxSsh(kx,f,N,params);
    So(:,t) = GmOmSsh(om,f,N,params);
    
    % from KE(om) to SSH(om)
    B=(freq1.'.^2-f.^2)./(freq1.'.^2+f.^2);
    SOobs(:,t)=A.*B.*ke;
    SKobs(:,t)=real(cg).'.*2.*SOobs(:,t);
    
    
%     hold on
%     loglog(freq,2*Sobs,'b')
%     loglog(om*86400/2/pi,So*2*pi/86400,'r')
%     hold off
    
    % now try to go from SSH(om) to SSH(kH)
    % SSH(kH)=cg*SSH(om)
    
%     figure;
%     hold on
%     loglog(kH*1e3,real(cg).'.*2.*Sobs,'b')
%     loglog(kx*1e3,Sx*1e-3,'r')
%     hold off
    
    
end

GMContinuum(i).kH{ii}{iii}    = kH;
GMContinuum(i).SOobs{ii}{iii} = SOobs;
GMContinuum(i).SKobs{ii}{iii} = SKobs;
GMContinuum(i).Sx{ii}{iii}    = Sx;
GMContinuum(i).So{ii}{iii}    = So;
GMContinuum(i).om{ii}{iii}    = om;
GMContinuum(i).kx{ii}{iii}    = kx;
GMContinuum(i).freq{ii}{iii}=freq;



                    end
                end
            end
        end
    end
    if rem(i,10)
        disp('save')
        save('../alb_mat/GMssh.mat','GMContinuum','-v7.3')
    end

end

  
%   
%   
%   % SSH frequency spectrum
%   subplot(3,1,1);
%   loglog(om*86400/2/pi,So*2*pi/86400,col);hold on;
%   %set(gca,'xlim',[f N],'xtick',10.^[-5:2]);
%   set(gca,'xlim',[0.8*f 1.2*N]*86400/2/pi,'xtick',10.^[0:2]);
%   ylabel('S_{\eta}(\omega) [m^2 (cpd)^{-1}]');
%   xlabel('\omega [cpd]');
%   grid on
%   ax = axis;
%   plot(f*[1 1]*86400/2/pi,[ax(3) ax(4)],'r--',...
%       N*[1 1]*86400/2/pi,[ax(3) ax(4)],'r--');
%   text(f*86400/2/pi*1.1,ax(3)*3,'f','color','r');
%   text(N*86400/2/pi*0.9,ax(3)*3,'N','color','r');
%   title('SSH Spectra');
%   
%   % SSH horizontal wavenumber spectrum
%   subplot(3,1,2);
%   loglog(kx*1e3,Sx*1e-3,col);hold on;
%   set(gca,'xlim',10.^[-3 2],'xtick',10.^[-3:2]);
%   ylabel('S_{\eta}(k_x) [m^2 cpkm^{-1}]');
%   xlabel('k_x [cpkm]');
%   grid on
%   xlin = [1:30]/10;
%   plot(xlin,xlin.^(-2)/xlin(1)^(-2)*10,'r--');
%   text(0.3,10,'-2 slope','color','r');
%   
%   % Horizontal wavenumber spectrum of surface slope
%   subplot(3,1,3);
%   loglog(kx*1e3,pdif(kx,Sx)*1e-3,col);hold on;
%   set(gca,'xlim',10.^[-3 2],'xtick',10.^[-3:2]);
%   ylabel('S_{\eta_x}(k_x) [cpkm^{-1}]');
%   xlabel('k_x [cpkm]');
%   grid on
% 
