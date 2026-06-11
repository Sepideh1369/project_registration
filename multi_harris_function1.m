# create multi_harris_function1


function [harris_function,gradient,angle]=multi_harris_function1(image,Mmax,ratio,sigma,d,perc,nbin,is_auto)

[M,N]=size(image);

%[M,N,P]=size(nonelinear_space);
harris_function=zeros(M,N,Mmax);
gradient=zeros(M,N,Mmax);
angle=zeros(M,N,Mmax);
Mask = image==0.001; SE = strel('square',5); Mask = imdilate(Mask,SE);Mask=1-Mask;
for i=1:Mmax
    scale=sigma*ratio^(i-1);
    radius=round(2*scale);
    j=-radius:1:radius;
    k=-radius:1:radius;
    [xarry,yarry]=meshgrid(j,k);
    if(strcmp(is_auto,'NO'))
        [k_percentile]=compute_k_percentile(xarry,yarry,perc,nbin);
    elseif(strcmp(is_auto,'YES'))
        [k_percentile]=compute_k_percentile_auto(xarry,yarry,perc);
    end

    %W=exp(-(abs(xarry)+abs(yarry))/scale);
   % k=1.2;
     W=1-exp(-3.315./((xarry.^2+yarry.^2).^4/k_percentile^8));
     W34=zeros(2*radius+1,2*radius+1);
     W12=zeros(2*radius+1,2*radius+1);
     W14=zeros(2*radius+1,2*radius+1);
     W23=zeros(2*radius+1,2*radius+1);
     W34(radius+2:2*radius+1,:)=W(radius+2:2*radius+1,:);
     W12(1:radius,:)=W(1:radius,:);
     W14(:,radius+2:2*radius+1)=W(:,radius+2:2*radius+1);
     W23(:,1:radius)=W(:,1:radius);
     M34=imfilter(image,W34,'replicate');
     M12=imfilter(image,W12,'replicate');
     M14=imfilter(image,W14,'replicate');
     M23=imfilter(image,W23,'replicate');

     Gx=log(M14./M23); Gx=Gx.*Mask;
     Gy=log(M34./M12); Gy=Gy.*Mask;

    Gx(find(imag(Gx)))=abs(Gx(find(imag(Gx))));
    Gy(find(imag(Gy)))=abs(Gy(find(imag(Gy))));
    Gx(~isfinite(Gx))=0;
    Gy(~isfinite(Gy))=0;
    temp_gradient=sqrt(Gx.^2+Gy.^2);
    temp_gradient=temp_gradient/max(temp_gradient(:));
    gradient(:,:,i)= temp_gradient;
    temp_angle=atan(Gy./Gx); temp_angle(isnan(temp_angle))=0;
    temp_angle=temp_angle/pi*180;
    temp_angle(temp_angle<0)=temp_angle(temp_angle<0)+180;
    angle(:,:,i)=temp_angle;
    Csh_11=scale^2*Gx.^2;
    Csh_12=scale^2*Gx.*Gy;
    Csh_22=scale^2*Gy.^2;
    gaussian_sigma=sqrt(2)*scale;
    width=round(3*gaussian_sigma);
    width_windows=2*width+1;
    W_gaussian=fspecial('gaussian',[width_windows width_windows],gaussian_sigma);
   [a,b]=meshgrid(1:width_windows,1:width_windows);
   index=find(((a-width-1)^2+(b-width-1)^2)>width^2);
   W_gaussian(index)=0;
   Csh_11=imfilter(Csh_11,W_gaussian,'replicate');
   Csh_12=imfilter(Csh_12,W_gaussian,'replicate');
   Csh_21=Csh_12;
   Csh_22=imfilter(Csh_22,W_gaussian,'replicate');
harris_function(:,:,i)=Csh_11.*Csh_22-Csh_21.*Csh_12-d*(Csh_11+Csh_22).^2;

%addpath('C:\Users\sepideh\Documents\MATLAB\Harris_Anisotropic scale space');

end

