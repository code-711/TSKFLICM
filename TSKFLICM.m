function [I_TSKFLICM,TSKFLICM_Vpe,TSKFLICM_SA,TSKFLICM_psnr]=TSKFLICM(data,v1,m,n,c,mc,e,ct)




[u,v1,Wij]=KWFLICM(data,v1,m,n,c,mc,e,ct);

ave=data;
t1=u;
G=zeros(m,n,c);
 p=zeros(m,n,c);



 
    for i=2:m-1
        for j=2:n-1
            tp1=0.0;
            for i1=-1:1
                for j1=-1:1
                   tp1=tp1+data(i+i1,j+j1);
                end
            end
             ave(i,j)=floor(tp1/9);  
             
        end
    end



   for i=1:m
        for j=1:n
            A(i,j)=sqrt(1+4*(a*ave(i,j)*ave(i,j)+b)^dd);
        end
   end


 
for k=1:c
    tp1=0.0;
    tp2=0.0;
    for i=1:m
        for j=1:n
            tp1=tp1+(u(i,j,k)^mc)*(t1(i,j,k)^nc)*(((a*data(i,j)*data(i,j)+b)^dd-2*(a*data(i,j)*v1(k)+b)^dd+(a*v1(k)*v1(k)+b)^dd)/A(i,j));
            tp2=tp2+(u(i,j,k)^mc)*(t1(i,j,k)^nc);
        end
    end
    gamma(k)=2*tp1/tp2;
end
  
 
  while et>0.0001 && t<1000 
    v=v1;

for k=1:c
    for i=2:m-1
        for j=2:n-1
            tp1=0.0;
            for i1=-1:1
                for j1=-1:1
                   if i1~=0 || j1~=0
                       tp1=tp1+exp(-(i1^2+j1^2)/8)*u(i+i1,j+j1,k);
                   end
                end
            end
           p(i,j,k)=tp1;
       end
    end
 end
 for i=1:m
        for j=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+(p(i,j,k))^f+0.0001;
            end
            for k=1:c
                pai(i,j,k)=((p(i,j,k))^f)/tp1;
            end
        end
 end




 for k=1:c
    for i=2:m-1
        for j=2:n-1
            tp1=0.0;
            for i1=-1:1
                for j1=-1:1
                   if i1~=0 || j1~=0
                    tp1=tp1+(1-pai(i+i1,j+j1,k))*Wij(i+i1,j+j1)*(1-u(i+i1,j+j1,k))^mc*(1-t1(i+i1,j+j1,k))^nc*((a*data(i+i1,j+j1)*data(i+i1,j+j1)+b)^dd-2*(a*data(i+i1,j+j1)*v1(k)+b)^dd+(a*v1(k)*v1(k)+b)^dd)/A(i+i1,j+j1);
                   end
                end
            end
            G(i,j,k)=tp1;
       end
    end
 end
 

 for k=1:c
        for  i=1:m
            for j=1:n
                d(i,j,k)=(1-pai(i,j,k))*(((a*data(i,j)*data(i,j)+b)^dd-2*(a*data(i,j)*v1(k)+b)^dd+(a*v1(k)*v1(k)+b)^dd))/A(i,j)+G(i,j,k)+0.0001;
            end
        end
 end


    for i=1:m
        for j=1:n
            tp1=0.0;
            for k=1:c
                tp1=tp1+((t1(i,j,k))^(nc-1)*d(i,j,k))^(-1/(mc-1));
            end
            for k=1:c
                u(i,j,k)=(((t1(i,j,k))^(nc-1)*d(i,j,k))^(-1/(mc-1)))/tp1;
            end
        end
    end
    
  
    for k=1:c
        for i=1:m
            for j=1:n
                 t1(i,j,k)=1/(1+(d(i,j,k)/((1-pai(i,j,k))*gamma(k)))^(1/(nc-1)));
            end
        end
    end
    
    for k=1:c
        tp1=0.0;
        tp2=0.0;
        for i=1:m
            for j=1:n
                tp1=tp1+((1-pai(i,j,k))*u(i,j,k)^mc*t1(i,j,k)^nc*(a*data(i,j)*v1(k)+b)^(dd-1)*data(i,j))/A(i,j);
                tp2=tp2+((1-pai(i,j,k))*u(i,j,k)^mc*t1(i,j,k)^nc*(a*v1(k)*v1(k)+b)^(dd-1))/A(i,j);
            end
        end
        v1(k)=tp1/tp2;            
    end


   temp=0.0;
   for k=1:c
         temp=temp+(v(k)-v1(k))^2;
   end
   if   temp<0.0001
        et=0.0001;
   end
t=t+1;
  end




I_TSKFLICM=zeros(m,n);
if c==2
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)
                I_TSKFLICM(i,j)=0;
            else
                I_TSKFLICM(i,j)=255;
            end
        end
    end
elseif c==3
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)
                I_TSKFLICM(i,j)=50;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)
                I_TSKFLICM(i,j)=100;
            else
                I_TSKFLICM(i,j)=255;
            end
        end
    end
elseif c==4
   for i=1:m
        for j=1:n
             if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)
                I_TSKFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)
                I_TSKFLICM(i,j)=60;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)
                I_TSKFLICM(i,j)=160;
            else
                I_TSKFLICM(i,j)=255;
            end
        end
    end
elseif c==5
    for i=1:m
        for j=1:n
            if u(i,j,1)>u(i,j,2)&& u(i,j,1)>u(i,j,3)&& u(i,j,1)>u(i,j,4)&&u(i,j,1)>u(i,j,5)
                I_TSKFLICM(i,j)=0;
            elseif u(i,j,2)>u(i,j,3)&& u(i,j,2)>u(i,j,1)&&u(i,j,2)>u(i,j,4)&&u(i,j,2)>u(i,j,5)
                I_TSKFLICM(i,j)=100;
            elseif u(i,j,3)>u(i,j,1)&& u(i,j,3)>u(i,j,2)&&u(i,j,3)>u(i,j,4)&&u(i,j,3)>u(i,j,5)
                I_TSKFLICM(i,j)=155;
            elseif u(i,j,4)>u(i,j,1)&& u(i,j,4)>u(i,j,2)&&u(i,j,4)>u(i,j,3)&&u(i,j,4)>u(i,j,5)
                I_TSKFLICM(i,j)=200;
            else
                I_TSKFLICM(i,j)=255;
            end
        end
    end
end
I=I_TSKFLICM;
[TSKFLICM_Vpe,TSKFLICM_SA,TSKFLICM_psnr]=huafenzhibiao(u,m,n,I,c);



