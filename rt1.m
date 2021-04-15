%Date: Feb. 16, 2021
%Author: Minah Oh
%Institution: James Madison University
%Email: ohmx@jmu.edu
%
% You can access the basis functions for the Order 1 Raviart Thomas Space (RT1) 
% in RTBasis. This is a three-dimensional array of size
% 8-by-8-by-NumberOfTriangles. 
%
% 1. Domain information is saved in meshgeometry.m
% 2. The corresponding "new_ele" files for mesh levels 1 through 9 are also provided here. 
% 3. "new_ele" files provide the three edge numbers for each triangle. This
%    determines the local ordering of the edges used in item 5 below.
%
% 4. Structure of the variable "RTBasis". RTBasis(:,:,k) provides the information 
%    needed to get the eight basis functions for the k-th triangle. 
%    For i=1, 2, ..., 8, 
%    basisX=@(x,y,i) RTBasis(1,i,k).*x+RTBasis(2,i,k).*y+RTBasis(3,i,k)+RTBasis(7,i,k).*x.*x+RTBasis(8,i,k).*x.*y;
%    basisY=@(x,y,i) RTBasis(4,i,k).*x+RTBasis(5,i,k).*y+RTBasis(6,i,k)+RTBasis(7,i,k).*x.*y+RTBasis(8,i,k).*y.*y;
%    where basisX and basisY are the X and Y coordinates of the i-th basis
%    function respectively. 
%
% 5. RTBasis(:,1:2,k) represents the two basis functions corresponding to the first edge of the k-th triangle.
%    RTBasis(:,3:4,k) represents the two basis functions corresponding to the second edge of the k-th triangle.
%    RTBasis(:,5:6,k) represents the two basis functions corresponding to the third edge of the k-th triangle.
%    RTBasis(:,7:8,k) represents the two basis functions corresponding to
%    the k-th triangle.
%
% 6. Global indexing of the basis functions: 
%    (1) 1~NumEdges: 
%        the first basis function corresponding to each edge.
%        (columns 1,3, and 5 of RTBasis(:,:,k))
%    (2) NumEdges+1~2*NumEdges: 
%        the second basis function corresponding to each edge
%        (columns 2,4, and 6 of RTBasis(:,:,k))
%    (3) 2*NumEdges+1~2*NumEdges+NumTriangles: 
%        the first basis function corresponding to each triangle.
%        (column 7 of RTBasis(:,:,k))
%    (4) 2*NumEdges+NumTriangles+1~2*NumEdges+2*NumTriangles: 
%        the second basis function corresponding to each triangle.
%        (column 8 of RTBasis(:,:,k))
%
% 7. Input mesh_level below and run the program. The output is the L2-error
%    of the toy problem used to check correctness of the code. 


%Input:
mesh_level=2;

exactx=@(x,y) 1./3.*x.^3.*y.^2-0.5.*x.^2.*y.^2;
exacty=@(x,y) -0.5.*x.^2.*y.^2+0.5.*x.*y.^2;

exact_div=@(x,y) x.^2.*y.^2 - x.^2.*y - x.*y.^2 + x.*y;

RHSx=@(x,y) exactx(x,y) - (y - 2.*x.*y + 2.*x.*y.^2 - y.^2);
RHSy=@(x,y) exacty(x,y) - (x - 2.*x.*y + 2.*x.^2.*y - x.^2);


model=createpde(1);
[a,b,c]=meshgeometry(1);
g=decsg(a,b,c);
geometryFromEdges(model,g);
[p,e,t]=initmesh(g,'hmax',inf);

for ii=1:mesh_level-1
   [p,e,t]=refinemesh(g,p,e,t,'regular'); 
   %pdemesh(p,e,t);
end

[~,N_node]=size(p);
node=p';
[~,N_ele]=size(t);
ele=t(1:3,1:N_ele);
ele=ele';
TR=triangulation(ele,node);
edge=edges(TR);
[N_edge,E2]=size(edge);
%load(['new_ele',num2str(mesh_level),'.mat']);
aa=load(['new_ele',num2str(mesh_level),'.mat']);
new_ele=cell2mat(struct2cell(aa));

I=zeros(N_ele*64,1);
J=zeros(N_ele*64,1);
S=zeros(N_ele*64,1);

i2=zeros(N_ele*8,1);
s2=zeros(N_ele*8,1);

RTBasis=zeros(8,8,N_ele);
abc=ones(8,8);
P1Basis=zeros(3,3,N_ele);

count=1;
count2=1;

for k=1:N_ele
   
    xx=[node(ele(k,1),1),node(ele(k,2),1),node(ele(k,3),1)];
    yy=[node(ele(k,1),2),node(ele(k,2),2),node(ele(k,3),2)];
    area=polyarea(xx,yy);
    
    %Check orientation of edges
    V=[1, node(ele(k,1),1), node(ele(k,1),2);...
       1, node(ele(k,2),1), node(ele(k,2),2);...
       1, node(ele(k,3),1), node(ele(k,3),2)];
    dV=det(V);
    if dV>0
       orientation=1; %Counterclockwise
    else 
       orientation=-1;
    end
    
    %Determine what type of triangle-upside down right triangle or standing
    %up right triangle
    if node(ele(k,1),2)==node(ele(k,2),2) && node(ele(k,3),2)>node(ele(k,1),2)
           Ttype=1;
           if node(ele(k,1),1)<node(ele(k,2),1)
              x1=node(ele(k,1),1);
              y1=node(ele(k,1),2);
              x2=node(ele(k,2),1);
              y2=node(ele(k,2),2);
              x3=node(ele(k,3),1);
              y3=node(ele(k,3),2);
           else
              x1=node(ele(k,2),1);
              y1=node(ele(k,2),2);
              x2=node(ele(k,1),1);
              y2=node(ele(k,1),2);
              x3=node(ele(k,3),1);
              y3=node(ele(k,3),2); 
           end
       
    elseif node(ele(k,1),2)==node(ele(k,3),2) && node(ele(k,2),2)>node(ele(k,1),2)
           Ttype=1;
             if node(ele(k,1),1)<node(ele(k,3),1)
              x1=node(ele(k,1),1);
              y1=node(ele(k,1),2);
              x2=node(ele(k,3),1);
              y2=node(ele(k,3),2);
              x3=node(ele(k,2),1);
              y3=node(ele(k,2),2);
           else
              x1=node(ele(k,3),1);
              y1=node(ele(k,3),2);
              x2=node(ele(k,1),1);
              y2=node(ele(k,1),2);
              x3=node(ele(k,2),1);
              y3=node(ele(k,2),2); 
             end
           
    elseif node(ele(k,2),2)==node(ele(k,3),2) && node(ele(k,1),2)>node(ele(k,2),2)
           Ttype=1;
             if node(ele(k,2),1)<node(ele(k,3),1)
              x1=node(ele(k,2),1);
              y1=node(ele(k,2),2);
              x2=node(ele(k,3),1);
              y2=node(ele(k,3),2);
              x3=node(ele(k,1),1);
              y3=node(ele(k,1),2);
           else
              x1=node(ele(k,3),1);
              y1=node(ele(k,3),2);
              x2=node(ele(k,2),1);
              y2=node(ele(k,2),2);
              x3=node(ele(k,1),1);
              y3=node(ele(k,1),2); 
           end
    else
        Ttype=2;
        for i=1:3
            k1=mod(i+1,3);
            if k1==0
               k1=3; 
            end
             k2=mod(i+2,3);
            if k2==0
               k2=3; 
            end
            
           if node(ele(k,i),2)<node(ele(k,k1),2)
              x1=node(ele(k,i),1);
              y1=node(ele(k,i),2);
              if node(ele(k,k1),1)<node(ele(k,k2),1)
                 x3=node(ele(k,k1),1);
                 y3=node(ele(k,k1),2);
                 x2=node(ele(k,k2),1);
                 y2=node(ele(k,k2),2);
              else 
                 x2=node(ele(k,k1),1);
                 y2=node(ele(k,k1),2);
                 x3=node(ele(k,k2),1);
                 y3=node(ele(k,k2),2); 
              end
           end
        end
    end
  
  for j=1:3
  xs=node(edge(new_ele(k,j),1),1);
  xe=node(edge(new_ele(k,j),2),1);
  ys=node(edge(new_ele(k,j),1),2);
  ye=node(edge(new_ele(k,j),2),2);
  
  if orientation<0
     dum=xs;
     xs=xe;
     xe=dum;
     dum=ys;
     ys=ye;
     ye=dum;
  end
  
  edge_length=sqrt((xs-xe)^2+(ys-ye)^2);
  %edge_length=1;
  
  if xs==xe 
      
      abc(2*j-1,1)=edge_length.*(xe/2 + xs/2);
      abc(2*j-1,2)=edge_length.*(ye/2 + ys/2);
      abc(2*j-1,3)=edge_length;
      abc(2*j-1,4)=0;
      abc(2*j-1,5)=0;
      abc(2*j-1,6)=0;
      abc(2*j-1,7)=edge_length.*(xe.^2./3 + (xe.*xs)./3 + xs.^2./3);
      abc(2*j-1,8)=edge_length.*((xe.*ye)./3 + (xe.*ys)./6 + (xs.*ye)./6 + (xs.*ys)./3);

      abc(2*j,1)=edge_length.*(xe./3 + xs./6);
      abc(2*j,2)=edge_length.*(ye./3 + ys./6);
      abc(2*j,3)=edge_length./2;
      abc(2*j,4)=0;
      abc(2*j,5)=0;
      abc(2*j,6)=0;
      abc(2*j,7)=edge_length.*(xe.^2./4 + (xe.*xs)./6 + xs.^2./12);
      abc(2*j,8)=edge_length.*((xe.*ye)./4 + (xe.*ys)./12 + (xs.*ye)./12 + (xs.*ys)./12);
      
     
   
  elseif ys==ye
      
     abc(2*j-1,1)=0;
     abc(2*j-1,2)=0;
     abc(2*j-1,3)=0;
     abc(2*j-1,4)=edge_length.*(xe./2 + xs./2);
     abc(2*j-1,5)=edge_length.*(ye./2 + ys./2);
     abc(2*j-1,6)=edge_length;
     abc(2*j-1,7)=edge_length.*((xe.*ye)./3 + (xe.*ys)./6 + (xs.*ye)./6 + (xs.*ys)./3);
     abc(2*j-1,8)=edge_length.*(ye.^2./3 + (ye.*ys)./3 + ys.^2./3);
     
     abc(2*j,1)=0;
     abc(2*j,2)=0;
     abc(2*j,3)=0;
     abc(2*j,4)= edge_length.*(xe./3 + xs./6);
     abc(2*j,5)=edge_length.*(ye./3 + ys./6);
     abc(2*j,6)=edge_length./2;
     abc(2*j,7)=edge_length.*((xe.*ye)./4 + (xe.*ys)./12 + (xs.*ye)./12 + (xs.*ys)./12);
     abc(2*j,8)=edge_length.*(ye.^2./4 + (ye.*ys)./6 + ys.^2./12);
     
  else
      
      abc(2*j-1,1)=edge_length.*((2.^(1./2).*xe)/4 + (2.^(1/2).*xs)./4);
      abc(2*j-1,2)=edge_length.*((2.^(1./2).*ye)./4 + (2.^(1/2).*ys)./4);
      abc(2*j-1,3)=(2.^(1./2).*edge_length)./2;
      abc(2*j-1,4)=edge_length.*((2.^(1./2).*xe)./4 + (2.^(1./2).*xs)./4);
      abc(2*j-1,5)=edge_length.*((2.^(1./2).*ye)./4 + (2.^(1./2).*ys)./4);
      abc(2*j-1,6)=(2.^(1./2).*edge_length)./2;
      abc(2*j-1,7)=edge_length.*((2.^(1./2).*xe.^2)./6 + (2.^(1./2).*xs.^2)./6 + (2.^(1./2).*xe.*xs)./6 + (2.^(1./2).*xe.*ye)./6 + (2.^(1./2).*xe.*ys)./12 + (2.^(1./2).*xs.*ye)./12 + (2.^(1./2).*xs.*ys)./6);
      abc(2*j-1,8)=edge_length*((2^(1/2)*ye^2)/6 + (2^(1/2)*ys^2)/6 + (2^(1/2)*xe*ye)/6 + (2^(1/2)*xe*ys)/12 + (2^(1/2)*xs*ye)/12 + (2^(1/2)*xs*ys)/6 + (2^(1/2)*ye*ys)/6);

      abc(2*j,1)=edge_length*((2^(1/2)*xe)/6 + (2^(1/2)*xs)/12);
      abc(2*j,2)=edge_length*((2^(1/2)*ye)/6 + (2^(1/2)*ys)/12);
      abc(2*j,3)=(2^(1/2)*edge_length)/4;
      abc(2*j,4)=edge_length*((2^(1/2)*xe)/6 + (2^(1/2)*xs)/12);
      abc(2*j,5)=edge_length*((2^(1/2)*ye)/6 + (2^(1/2)*ys)/12);
      abc(2*j,6)=(2^(1/2)*edge_length)/4;
      abc(2*j,7)=edge_length*((2^(1/2)*xe^2)/8 + (2^(1/2)*xs^2)/24 + (2^(1/2)*xe*xs)/12 + (2^(1/2)*xe*ye)/8 + (2^(1/2)*xe*ys)/24 + (2^(1/2)*xs*ye)/24 + (2^(1/2)*xs*ys)/24);
      abc(2*j,8)=edge_length*((2^(1/2)*ye^2)/8 + (2^(1/2)*ys^2)/24 + (2^(1/2)*xe*ye)/8 + (2^(1/2)*xe*ys)/24 + (2^(1/2)*xs*ye)/24 + (2^(1/2)*xs*ys)/24 + (2^(1/2)*ye*ys)/12);

  end
  
  end
  
  if Ttype==1
  abc(7,1)=x1^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - x2^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - (x1^3*(y2 - y3))/(3*(x2 - x3)) + (x2^3*(y2 - y3))/(3*(x2 - x3));

  abc(7,2)=(x2^3*(y2 - y3)^2)/(6*(x2 - x3)^2) - (x1^3*(y2 - y3)^2)/(6*(x2 - x3)^2) + (x1*(x2*y1 - x3*y1 + x2*y3 - x3*y2)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2) - (x2*(x2*y1 - x3*y1 + x2*y3 - x3*y2)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2) - (x1^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2) + (x2^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2);
 
  abc(7,3)=(x2^2*(y2 - y3))/(2*(x2 - x3)) - (x1^2*(y2 - y3))/(2*(x2 - x3)) + (x1*(2*x2 - 2*x3)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2) - (x2*(2*x2 - 2*x3)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2);

  abc(7,4)=0;
  abc(7,5)=0;
  abc(7,6)=0;
  abc(7,7)=x1^3*(y1/3 - (x2*y3 - x3*y2)/(3*(x2 - x3))) - x2^3*(y1/3 - (x2*y3 - x3*y2)/(3*(x2 - x3))) - (x1^4*(2*x2 - 2*x3)*(y2 - y3))/(8*(x2 - x3)^2) + (x2^4*(2*x2 - 2*x3)*(y2 - y3))/(8*(x2 - x3)^2);

  abc(7,8)=x2^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - x1^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - (x1^4*(y2 - y3)^2)/(8*(x2 - x3)^2) + (x2^4*(y2 - y3)^2)/(8*(x2 - x3)^2) - (x1^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2) + (x2^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2);
 


  abc(8,1)=0;
  abc(8,2)=0;
  abc(8,3)=0;
  abc(8,4)=x1^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - x2^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - (x1^3*(y2 - y3))/(3*(x2 - x3)) + (x2^3*(y2 - y3))/(3*(x2 - x3));

  abc(8,5)=x2*((x2*y3 - x3*y2)^2/(2*(x2 - x3)^2) - y1^2/2) - x1*((x2*y3 - x3*y2)^2/(2*(x2 - x3)^2) - y1^2/2) - (x1^3*(y2 - y3)^2)/(6*(x2 - x3)^2) + (x2^3*(y2 - y3)^2)/(6*(x2 - x3)^2) - (x1^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2) + (x2^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2);

  abc(8,6)=x1*(y1 - (x2*y3 - x3*y2)/(x2 - x3)) - x2*(y1 - (x2*y3 - x3*y2)/(x2 - x3)) - (x1^2*(y2 - y3))/(2*(x2 - x3)) + (x2^2*(y2 - y3))/(2*(x2 - x3));

  abc(8,7)=x2^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - x1^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - (x1^4*(3*x2 - 3*x3)*(y2 - y3)^2)/(24*(x2 - x3)^3) + (x2^4*(3*x2 - 3*x3)*(y2 - y3)^2)/(24*(x2 - x3)^3) - (x1^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2) + (x2^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2);
 
  abc(8,8)=x2*((x2*y3 - x3*y2)^3/(3*(x2 - x3)^3) - y1^3/3) - x1*((x2*y3 - x3*y2)^3/(3*(x2 - x3)^3) - y1^3/3) - (x1^2*(x2*y3 - x3*y2)^2*(y2 - y3))/(2*(x2 - x3)^3) - (x1^3*(x2*y3 - x3*y2)*(y2 - y3)^2)/(3*(x2 - x3)^3) + (x2^2*(x2*y3 - x3*y2)^2*(y2 - y3))/(2*(x2 - x3)^3) + (x2^3*(x2*y3 - x3*y2)*(y2 - y3)^2)/(3*(x2 - x3)^3) - (x1^4*(2*y2 - 2*y3)*(y2 - y3)^2)/(24*(x2 - x3)^3) + (x2^4*(2*y2 - 2*y3)*(y2 - y3)^2)/(24*(x2 - x3)^3);
 

  
  else
  abc(7,1)=((x2 - x3)^2*(y1 - y3)*(4*x3^2 + 8*x2*x3))/(24*(x1 - x3)^2) - (x1*(8*x2 + 4*x3)*(x2 - x3)^2*(y1 - y3))/(24*(x1 - x3)^2);
  abc(7,2)=((x2 - x3)^2*(y1 - y3)*(4*x3*y1 - 4*x2*y1 + 4*x2*y3 + 8*x3*y3))/(24*(x1 - x3)^2) - (x1*y3*(x2 - x3)^2*(y1 - y3))/(2*(x1 - x3)^2);
  abc(7,3)=(x3*(x2 - x3)^2*(y1 - y3))/(2*(x1 - x3)^2) - (x1*(x2 - x3)^2*(y1 - y3))/(2*(x1 - x3)^2);
  abc(7,4)=0;
  abc(7,5)=0;
  abc(7,6)=0;
  abc(7,7)=((x2 - x3)^2*(y1 - y3)*(6*x2^2*x3 + 4*x2*x3^2 + 2*x3^3))/(24*(x1 - x3)^2) - (x1*(x2 - x3)^2*(y1 - y3)*(6*x2^2 + 4*x2*x3 + 2*x3^2))/(24*(x1 - x3)^2);
  abc(7,8)=((x2 - x3)^2*(y1 - y3)*(x3^2*y1 - 3*x2^2*y1 + 3*x2^2*y3 + 3*x3^2*y3 + 2*x2*x3*y1 + 6*x2*x3*y3))/(24*(x1 - x3)^2) - (x1*(8*x2*y3 + 4*x3*y3)*(x2 - x3)^2*(y1 - y3))/(24*(x1 - x3)^2);    
  
  abc(8,1)=0;
  abc(8,2)=0;
  abc(8,3)=0;
  abc(8,4)=-((x2 - x3)^2*(y1 - y3)*(4*x1^2*x3 + 8*x2*x1^2 - 8*x1*x3^2 - 16*x2*x1*x3 + 4*x3^3 + 8*x2*x3^2))/(24*(x1 - x3)^3);
  abc(8,5)= -((x2 - x3)^2*(y1 - y3)*(12*x1^2*y3 + 4*x3^2*y1 + 8*x3^2*y3 + 4*x1*x2*y1 - 4*x1*x3*y1 - 4*x1*x2*y3 - 4*x2*x3*y1 - 20*x1*x3*y3 + 4*x2*x3*y3))/(24*(x1 - x3)^3);
  abc(8,6)=-((x2 - x3)^2*(y1 - y3)*(12*x1^2 - 24*x1*x3 + 12*x3^2))/(24*(x1 - x3)^3);
  abc(8,7)=-((x2 - x3)^2*(y1 - y3)*(x3^3*y1 + 3*x3^3*y3 + 3*x1*x2^2*y1 - x1*x3^2*y1 - 3*x1*x2^2*y3 + 2*x2*x3^2*y1 + 8*x1^2*x2*y3 - 3*x2^2*x3*y1 - 7*x1*x3^2*y3 + 4*x1^2*x3*y3 + 6*x2*x3^2*y3 + 3*x2^2*x3*y3 - 2*x1*x2*x3*y1 - 14*x1*x2*x3*y3))/(24*(x1 - x3)^3);
  abc(8,8)=-((x2 - x3)^2*(y1 - y3)*(12*x1^2*y3^2 + 8*x1*x2*y1*y3 - 8*x1*x2*y3^2 - 8*x1*x3*y1*y3 - 16*x1*x3*y3^2 + 2*x2^2*y1^2 - 4*x2^2*y1*y3 + 2*x2^2*y3^2 - 4*x2*x3*y1^2 + 4*x2*x3*y3^2 + 2*x3^2*y1^2 + 4*x3^2*y1*y3 + 6*x3^2*y3^2))/(24*(x1 - x3)^3);
  end
  
%abc(7,:)=abc(7,:)./area;
%abc(8,:)=abc(8,:)./area;
  
 RTBasis(:,:,k)=abc\eye(8);
 

 %{
if Ttype==2
   for j=1:3
          xs=node(edge(new_ele(k,j),1),1);
          xe=node(edge(new_ele(k,j),2),1);
          ys=node(edge(new_ele(k,j),1),2);
          ye=node(edge(new_ele(k,j),2),2);
       if xs~=xe 
          if ys~=ye
          RTBasis(:,2*j-1,k) =-RTBasis(:,2*j-1,k);
          RTBasis(:,2*j,k)=-RTBasis(:,2*j,k);
          end
       end
   end
   
else
    for j=1:3
          xs=node(edge(new_ele(k,j),1),1);
          xe=node(edge(new_ele(k,j),2),1);
          ys=node(edge(new_ele(k,j),1),2);
          ye=node(edge(new_ele(k,j),2),2);
       if xs==xe
          RTBasis(:,2*j-1,k) =-RTBasis(:,2*j-1,k);
          RTBasis(:,2*j,k)=-RTBasis(:,2*j,k);
       end
       
       if ys==ye
          RTBasis(:,2*j-1,k) =-RTBasis(:,2*j-1,k);
          RTBasis(:,2*j,k)=-RTBasis(:,2*j,k);
       end
   end
    
end
%}

 
%RTBasis(:,:,k)
basisX=@(x,y,i) RTBasis(1,i,k).*x+RTBasis(2,i,k).*y+RTBasis(3,i,k)+RTBasis(7,i,k).*x.*x+RTBasis(8,i,k).*x.*y;
basisY=@(x,y,i) RTBasis(4,i,k).*x+RTBasis(5,i,k).*y+RTBasis(6,i,k)+RTBasis(7,i,k).*x.*y+RTBasis(8,i,k).*y.*y;
divBasis=@(x,y,i) RTBasis(1,i,k)+2.*RTBasis(7,i,k).*x+RTBasis(8,i,k).*y...
                 +RTBasis(5,i,k)+2.*RTBasis(8,i,k).*y+RTBasis(7,i,k).*x;

 [X,Y,Wx,Wy]=triquad(8,[node(ele(k,1),1),node(ele(k,1),2); node(ele(k,2),1),node(ele(k,2),2); node(ele(k,3),1),node(ele(k,3),2)]);

 %{
 g=@(x,y) basisX(x,y,7);
 check=Wx'*feval(g,X,Y)*Wy
 
 g=@(x,y) basisY(x,y,7);
 check=Wx'*feval(g,X,Y)*Wy
 
 g2=@(x,y) basisX(x,y,8);
 check=Wx'*feval(g2,X,Y)*Wy
 
 g2=@(x,y) basisY(X,Y,8);
 check=Wx'*feval(g2,X,Y)*Wy
%}
 
for s=1:8
    for t=s:8
       f=@(x,y) basisX(x,y,s).*basisX(x,y,t)+basisY(x,y,s).*basisY(x,y,t)+divBasis(x,y,s).*divBasis(x,y,t);
       %local(s,t,k)= Wx'*feval(f,X,Y)*Wy;
       S(count)=Wx'*feval(f,X,Y)*Wy;
       
       if s==1
          I(count)=new_ele(k,1);
       elseif s==2
          I(count)=N_edge+new_ele(k,1); 
       elseif s==3
          I(count)=new_ele(k,2);
       elseif s==4
          I(count)=N_edge+new_ele(k,2);
       elseif s==5
          I(count)=new_ele(k,3);
       elseif s==6
          I(count)=N_edge+new_ele(k,3);
       elseif s==7
          I(count)=2*N_edge+k;
       else 
          I(count)=2*N_edge+N_ele+k; 
       end
       
       
      if t==1
          J(count)=new_ele(k,1);
       elseif t==2
          J(count)=N_edge+new_ele(k,1); 
       elseif t==3
          J(count)=new_ele(k,2);
       elseif t==4
          J(count)=N_edge+new_ele(k,2);
       elseif t==5
          J(count)=new_ele(k,3);
       elseif t==6
          J(count)=N_edge+new_ele(k,3);
       elseif t==7
          J(count)=2*N_edge+k;
       else 
          J(count)=2*N_edge+N_ele+k; 
       end
    
    
    count=count+1;
    
    if s<t
       S(count)=S(count-1);
       I(count)=J(count-1);
       J(count)=I(count-1);
       count=count+1;
    end
    
    end

    f=@(x,y) RHSx(x,y).*basisX(x,y,s) + RHSy(x,y).*basisY(x,y,s);
    s2(count2)=Wx'*feval(f,X,Y)*Wy;
      if s==1
          i2(count2)=new_ele(k,1);
       elseif s==2
          i2(count2)=N_edge+new_ele(k,1); 
       elseif s==3
          i2(count2)=new_ele(k,2);
       elseif s==4
          i2(count2)=N_edge+new_ele(k,2);
       elseif s==5
          i2(count2)=new_ele(k,3);
       elseif s==6
          i2(count2)=N_edge+new_ele(k,3);
       elseif s==7
          i2(count2)=2*N_edge+k;
       else 
          i2(count2)=2*N_edge+N_ele+k; 
       end
       count2=count2+1;
end

end

matrix_size=N_edge+N_edge+N_ele+N_ele;
A=sparse(I,J,S,matrix_size,matrix_size);
b=sparse(i2,1,s2,matrix_size,1);

c=A\b;

%Computing L^2_r-error
error=0;

for k=1:N_ele
    
    [X,Y,Wx,Wy]=triquad(8,[node(ele(k,1),1),node(ele(k,1),2); node(ele(k,2),1),node(ele(k,2),2); node(ele(k,3),1),node(ele(k,3),2)]);
   
    basisX=@(x,y,i) RTBasis(1,i,k).*x+RTBasis(2,i,k).*y+RTBasis(3,i,k)+RTBasis(7,i,k).*x.*x+RTBasis(8,i,k).*x.*y;
    basisY=@(x,y,i) RTBasis(4,i,k).*x+RTBasis(5,i,k).*y+RTBasis(6,i,k)+RTBasis(7,i,k).*x.*y+RTBasis(8,i,k).*y.*y;

   uhx=@(x,y) c(new_ele(k,1)).*basisX(x,y,1)+c(N_edge+new_ele(k,1)).*basisX(x,y,2)...
             +c(new_ele(k,2)).*basisX(x,y,3)+c(N_edge+new_ele(k,2)).*basisX(x,y,4)...
             +c(new_ele(k,3)).*basisX(x,y,5)+c(N_edge+new_ele(k,3)).*basisX(x,y,6)...
               +c(2*N_edge+k).*basisX(x,y,7)...
               +c(2*N_edge+N_ele+k).*basisX(x,y,8);
   uhy=@(x,y) c(new_ele(k,1)).*basisY(x,y,1)+c(N_edge+new_ele(k,1)).*basisY(x,y,2)...
             +c(new_ele(k,2)).*basisY(x,y,3)+c(N_edge+new_ele(k,2)).*basisY(x,y,4)...
             +c(new_ele(k,3)).*basisY(x,y,5)+c(N_edge+new_ele(k,3)).*basisY(x,y,6)...
               +c(2*N_edge+k).*basisY(x,y,7)...
               +c(2*N_edge+N_ele+k).*basisY(x,y,8);
   f=@(x,y) (exactx(x,y)-uhx(x,y)).^2+(exacty(x,y)-uhy(x,y)).^2;    
   error=error+Wx'*feval(f,X,Y)*Wy;
end

error=sqrt(error);

%fprintf("L2 Error mesh level %d: %e",mesh_level,error);
%fprintf('\n');
