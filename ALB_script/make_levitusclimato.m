% from frequency spectrum to vertical wavenumber spectrum 
% sw_vmode /ARNAUD/SCRIPPS/IWISE/Matlablocal/MHA_MatlabLocal
% 
%
% 
libdir='~/ARNAUD/TOOLBOXES/';
addpath([libdir 'matlab_archive/ToolBox-MMP-CURRENT/trunk/seawater_ver3_2/'])
addpath(sprintf('%smatlab_archive/ToolBox-MMP-CURRENT/trunk/%s/',libdir,'Oceanography_toolbox'))

ncid=netcdf.open('/Users/aleboyer/ARNAUD/LEVITUS/data.nc');



[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncid);
%[varname,xtype,dimids,natts] = netcdf.inqVar(ncid,6)

time=netcdf.getVar(ncid,0)*30;
Y=netcdf.getVar(ncid,1);
Z=netcdf.getVar(ncid,2);
X=netcdf.getVar(ncid,3);
sal=netcdf.getVar(ncid,4);
temp=netcdf.getVar(ncid,5);

temp(temp<0)=NaN;
sal(sal<0)=NaN;

[LY,LZ,LX,LT]=size(temp);

[XX,YY]=meshgrid(X,Y);
ZZ=zeros(LZ,LY,LX);
for x=1:LX
    for y=1:LY
        ZZ(:,y,x)=pressure(Z,Y(y));
    end
end
sal1=permute(sal,[2 4 1 3]);
temp1=permute(temp,[2 4 1 3]);

sal1=reshape(sal1,[24 4 360*180]);
temp1=reshape(temp1,[24 4 360*180]);
ZZ=reshape(ZZ,[24 360*180]);
n2=zeros([23 4 360*180])*NaN;
for t=1:LT
    ta=squeeze(temp1(:,t,:));sa=squeeze(sal1(:,t,:));
    [n2(:,t,:),q,p_ave] = sw_bfrq(sa,ta,ZZ);
end

n2=reshape(n2,[23,4,180,360]);
p_ave=reshape(p_ave,[23,180,360]);


%% function to smmoth transition between season 
time1=[1; time; 365];
n2_1=zeros(23,4+2,180,360);
n2_1(:,2:end-1,:,:)=n2;
n2_1(:,1,:,:)=(n2(:,1,:,:)-n2(:,end,:,:))/(time(1)+(365-time(4)))*(365-time(4))+n2(:,end,:,:);
n2_1(:,end,:,:)=(n2(:,1,:,:)-n2(:,end,:,:))/(time(1)+(365-time(4)))*(364-time(4))+n2(:,end,:,:);


ydaytime=1:365;
mask=isnan(squeeze(n2_1(:,1,:,:)));
n2_2=zeros(23,365,180,360);
n2_1(isnan(n2_1))=0;
for y=1:LY
    for x=1:LX
        for z=1:LZ-1
            n2_2(z,:,y,x)=interp1(time1,squeeze(n2_1(z,:,y,x)),ydaytime,'pship');
        end
    end
end

climatology.n2 = n2_2;
climatology.yday = ydaytime;
climatology.lat  = Y;
climatology.lon  = X;
climatology.z    = p_ave;


save('levitus_N2climatology.mat','climatology','-v7.3')




