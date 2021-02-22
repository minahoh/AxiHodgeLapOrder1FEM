%Date: Feb. 17, 2021
%Author: Minah Oh
%Institution: James Madison University
%Email: ohmx@jmu.edu
%
% You can access the basis functions for the Order 1 Nedelec Space (ND1) 
% in NDBasis. This is a three-dimensional array of size
% 8-by-8-by-NumberOfTriangles. 
%
% 1. Domain information is saved in meshgeometry.m
% 2. The corresponding "new_ele" files for mesh levels 1 through 9 are also provided here. 
% 3. "new_ele" files provide the three edge numbers for each triangle. This
%    determines the local ordering of the edges used in item 5 below.
%
% 4. Structure of the variable "NDBasis". NDBasis(:,:,k) provides the information 
%    needed to get the eight basis functions for the k-th triangle. 
%    For i=1, 2, ..., 8, 
%    basisX=@(x,y,i) NDBasis(1,i,k).*x+NDBasis(2,i,k).*y+NDBasis(3,i,k)-NDBasis(7,i,k).*x.*y-NDBasis(8,i,k).*y.*y;
%    basisY=@(x,y,i) NDBasis(4,i,k).*x+NDBasis(5,i,k).*y+NDBasis(6,i,k)+NDBasis(7,i,k).*x.*x+NDBasis(8,i,k).*x.*y;
%    where basisX and basisY are the X and Y coordinates of the i-th basis
%    function respectively. 
%
% 5. NDBasis(:,1:2,k) represents the two basis functions corresponding to the first edge of the k-th triangle.
%    NDBasis(:,3:4,k) represents the two basis functions corresponding to the second edge of the k-th triangle.
%    NDBasis(:,5:6,k) represents the two basis functions corresponding to the third edge of the k-th triangle.
%    NDBasis(:,7:8,k) represents the two basis functions corresponding to
%    the k-th triangle.
%
% 6. Global indexing of the basis functions: 
%    (1) 1~NumEdges: 
%        the first basis function corresponding to each edge.
%        (columns 1,3, and 5 of NDBasis(:,:,k))
%    (2) NumEdges+1~2*NumEdges: 
%        the second basis function corresponding to each edge
%        (columns 2,4, and 6 of NDBasis(:,:,k))
%    (3) 2*NumEdges+1~2*NumEdges+NumTriangles: 
%        the first basis function corresponding to each triangle.
%        (column 7 of NDBasis(:,:,k))
%    (4) 2*NumEdges+NumTriangles+1~2*NumEdges+2*NumTriangles: 
%        the second basis function corresponding to each triangle.
%        (column 8 of NDBasis(:,:,k))
%
% 7. Input mesh_level below and run the program. The output is the L2-error
%    of the toy problem used to check correctness of the code. 


%Input:
mesh_level=7;


exactx=@(x,y) 10; 
exacty=@(x,y) -1/pi.*cos(pi.*x).*sin(pi.*y); 


RHSx=@(x,y) exactx(x,y)+pi.*sin(pi.*x).*cos(pi.*y);
RHSy=@(x,y) exacty(x,y)-pi.*cos(pi.*x).*sin(pi.*y);

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
load(['new_ele',num2str(mesh_level),'.mat']);

I=zeros(N_ele*64,1);
J=zeros(N_ele*64,1);
S=zeros(N_ele*64,1);

i2=zeros(N_ele*8,1);
s2=zeros(N_ele*8,1);

NDBasis=zeros(8,8,N_ele);
abc=ones(8,8);
P1Basis=zeros(3,3,N_ele);

count=1;
count2=1;

for k=1:N_ele
   
    %{
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
    %}
    
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
  
      
      abc(2*j-1,1)=xe^2/2 - xs^2/2;
      abc(2*j-1,2)=(xe*ye)/2 + (xe*ys)/2 - (xs*ye)/2 - (xs*ys)/2;
      abc(2*j-1,3)=xe - xs;
      abc(2*j-1,4)=(xe*ye)/2 - (xe*ys)/2 + (xs*ye)/2 - (xs*ys)/2;
      abc(2*j-1,5)=ye^2/2 - ys^2/2;
      abc(2*j-1,6)=ye - ys;
      abc(2*j-1,7)=(xs^2*ye)/2 - (xe^2*ys)/2 + (xe*xs*ye)/2 - (xe*xs*ys)/2;
      abc(2*j-1,8)=(xs*ye^2)/2 - (xe*ys^2)/2 - (xe*ye*ys)/2 + (xs*ye*ys)/2;

      abc(2*j,1)=xe^2/3 - (xe*xs)/6 - xs^2/6;
      abc(2*j,2)=(xe*ye)/3 + (xe*ys)/6 - (xs*ye)/3 - (xs*ys)/6;
      abc(2*j,3)=xe/2 - xs/2;
      abc(2*j,4)=(xe*ye)/3 - (xe*ys)/3 + (xs*ye)/6 - (xs*ys)/6;
      abc(2*j,5)=ye^2/3 - (ye*ys)/6 - ys^2/6;
      abc(2*j,6)=ye/2 - ys/2;
      abc(2*j,7)=(xs^2*ye)/6 - (xe^2*ys)/3 + (xe*xs*ye)/3 - (xe*xs*ys)/6;
      abc(2*j,8)=(xs*ye^2)/3 - (xe*ys^2)/6 - (xe*ye*ys)/3 + (xs*ye*ys)/6;
      
  end
  
  
  
  if Ttype==1
  
  abc(7,1)=x1^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - x2^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - (x1^3*(y2 - y3))/(3*(x2 - x3)) + (x2^3*(y2 - y3))/(3*(x2 - x3));

  abc(7,2)=x2*((x2*y3 - x3*y2)^2/(2*(x2 - x3)^2) - y1^2/2) - x1*((x2*y3 - x3*y2)^2/(2*(x2 - x3)^2) - y1^2/2) - (x1^3*(y2 - y3)^2)/(6*(x2 - x3)^2) + (x2^3*(y2 - y3)^2)/(6*(x2 - x3)^2) - (x1^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2) + (x2^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2);

  abc(7,3)=x1*(y1 - (x2*y3 - x3*y2)/(x2 - x3)) - x2*(y1 - (x2*y3 - x3*y2)/(x2 - x3)) - (x1^2*(y2 - y3))/(2*(x2 - x3)) + (x2^2*(y2 - y3))/(2*(x2 - x3));

  abc(7,4)=0;
  abc(7,5)=0;
  abc(7,6)=0;
  abc(7,7)=x1^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - x2^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) + (x1^4*(3*x2 - 3*x3)*(y2 - y3)^2)/(24*(x2 - x3)^3) - (x2^4*(3*x2 - 3*x3)*(y2 - y3)^2)/(24*(x2 - x3)^3) + (x1^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2) - (x2^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2);

  abc(7,8)=x1*((x2*y3 - x3*y2)^3/(3*(x2 - x3)^3) - y1^3/3) - x2*((x2*y3 - x3*y2)^3/(3*(x2 - x3)^3) - y1^3/3) + (x1^2*(x2*y3 - x3*y2)^2*(y2 - y3))/(2*(x2 - x3)^3) + (x1^3*(x2*y3 - x3*y2)*(y2 - y3)^2)/(3*(x2 - x3)^3) - (x2^2*(x2*y3 - x3*y2)^2*(y2 - y3))/(2*(x2 - x3)^3) - (x2^3*(x2*y3 - x3*y2)*(y2 - y3)^2)/(3*(x2 - x3)^3) + (x1^4*(2*y2 - 2*y3)*(y2 - y3)^2)/(24*(x2 - x3)^3) - (x2^4*(2*y2 - 2*y3)*(y2 - y3)^2)/(24*(x2 - x3)^3);
 

  abc(8,1)=0;
  abc(8,2)=0;
  abc(8,3)=0;
  abc(8,4)=x1^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - x2^2*(y1/2 - (x2*y3 - x3*y2)/(2*(x2 - x3))) - (x1^3*(y2 - y3))/(3*(x2 - x3)) + (x2^3*(y2 - y3))/(3*(x2 - x3));

  abc(8,5)=(x2^3*(y2 - y3)^2)/(6*(x2 - x3)^2) - (x1^3*(y2 - y3)^2)/(6*(x2 - x3)^2) + (x1*(x2*y1 - x3*y1 + x2*y3 - x3*y2)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2) - (x2*(x2*y1 - x3*y1 + x2*y3 - x3*y2)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2) - (x1^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2) + (x2^2*(x2*y3 - x3*y2)*(y2 - y3))/(2*(x2 - x3)^2);

  abc(8,6)=(x2^2*(y2 - y3))/(2*(x2 - x3)) - (x1^2*(y2 - y3))/(2*(x2 - x3)) + (x1*(2*x2 - 2*x3)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2) - (x2*(2*x2 - 2*x3)*(x2*y1 - x3*y1 - x2*y3 + x3*y2))/(2*(x2 - x3)^2);
 
  abc(8,7)=x1^3*(y1/3 - (x2*y3 - x3*y2)/(3*(x2 - x3))) - x2^3*(y1/3 - (x2*y3 - x3*y2)/(3*(x2 - x3))) - (x1^4*(2*x2 - 2*x3)*(y2 - y3))/(8*(x2 - x3)^2) + (x2^4*(2*x2 - 2*x3)*(y2 - y3))/(8*(x2 - x3)^2);

  abc(8,8)=x2^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - x1^2*((x2*y3 - x3*y2)^2/(4*(x2 - x3)^2) - y1^2/4) - (x1^4*(y2 - y3)^2)/(8*(x2 - x3)^2) + (x2^4*(y2 - y3)^2)/(8*(x2 - x3)^2) - (x1^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2) + (x2^3*(x2*y3 - x3*y2)*(y2 - y3))/(3*(x2 - x3)^2);


  
  else
  abc(7,1)=-((x2 - x3)^2*(y1 - y3)*(4*x1^2*x3 + 8*x2*x1^2 - 8*x1*x3^2 - 16*x2*x1*x3 + 4*x3^3 + 8*x2*x3^2))/(24*(x1 - x3)^3);
  abc(7,2)=-((x2 - x3)^2*(y1 - y3)*(12*x1^2*y3 + 4*x3^2*y1 + 8*x3^2*y3 + 4*x1*x2*y1 - 4*x1*x3*y1 - 4*x1*x2*y3 - 4*x2*x3*y1 - 20*x1*x3*y3 + 4*x2*x3*y3))/(24*(x1 - x3)^3); 
  abc(7,3)=-((x2 - x3)^2*(y1 - y3)*(12*x1^2 - 24*x1*x3 + 12*x3^2))/(24*(x1 - x3)^3);
  abc(7,4)=0;
  abc(7,5)=0;
  abc(7,6)=0;
  abc(7,7)=((x2 - x3)^2*(y1 - y3)*(x3^3*y1 + 3*x3^3*y3 + 3*x1*x2^2*y1 - x1*x3^2*y1 - 3*x1*x2^2*y3 + 2*x2*x3^2*y1 + 8*x1^2*x2*y3 - 3*x2^2*x3*y1 - 7*x1*x3^2*y3 + 4*x1^2*x3*y3 + 6*x2*x3^2*y3 + 3*x2^2*x3*y3 - 2*x1*x2*x3*y1 - 14*x1*x2*x3*y3))/(24*(x1 - x3)^3);
  abc(7,8)=((x2 - x3)^2*(y1 - y3)*(12*x1^2*y3^2 + 8*x1*x2*y1*y3 - 8*x1*x2*y3^2 - 8*x1*x3*y1*y3 - 16*x1*x3*y3^2 + 2*x2^2*y1^2 - 4*x2^2*y1*y3 + 2*x2^2*y3^2 - 4*x2*x3*y1^2 + 4*x2*x3*y3^2 + 2*x3^2*y1^2 + 4*x3^2*y1*y3 + 6*x3^2*y3^2))/(24*(x1 - x3)^3);
  
  abc(8,1)=0;
  abc(8,2)=0;  
  abc(8,3)=0;  
  abc(8,4)=((x2 - x3)^2*(y1 - y3)*(4*x3^2 + 8*x2*x3))/(24*(x1 - x3)^2) - (x1*(8*x2 + 4*x3)*(x2 - x3)^2*(y1 - y3))/(24*(x1 - x3)^2);  
  abc(8,5)=((x2 - x3)^2*(y1 - y3)*(4*x3*y1 - 4*x2*y1 + 4*x2*y3 + 8*x3*y3))/(24*(x1 - x3)^2) - (x1*y3*(x2 - x3)^2*(y1 - y3))/(2*(x1 - x3)^2);  
  abc(8,6)=(x3*(x2 - x3)^2*(y1 - y3))/(2*(x1 - x3)^2) - (x1*(x2 - x3)^2*(y1 - y3))/(2*(x1 - x3)^2);  
  abc(8,7)=((x2 - x3)^2*(y1 - y3)*(6*x2^2*x3 + 4*x2*x3^2 + 2*x3^3))/(24*(x1 - x3)^2) - (x1*(x2 - x3)^2*(y1 - y3)*(6*x2^2 + 4*x2*x3 + 2*x3^2))/(24*(x1 - x3)^2);  
  abc(8,8)=((x2 - x3)^2*(y1 - y3)*(x3^2*y1 - 3*x2^2*y1 + 3*x2^2*y3 + 3*x3^2*y3 + 2*x2*x3*y1 + 6*x2*x3*y3))/(24*(x1 - x3)^2) - (x1*(8*x2*y3 + 4*x3*y3)*(x2 - x3)^2*(y1 - y3))/(24*(x1 - x3)^2); 
  
  end
  
%abc(7,:)=abc(7,:)./area;
%abc(8,:)=abc(8,:)./area;
  
 NDBasis(:,:,k)=abc\eye(8);
 

%NDBasis(:,:,k)
basisX=@(x,y,i) NDBasis(1,i,k).*x+NDBasis(2,i,k).*y+NDBasis(3,i,k)-NDBasis(7,i,k).*x.*y-NDBasis(8,i,k).*y.*y; 
basisY=@(x,y,i) NDBasis(4,i,k).*x+NDBasis(5,i,k).*y+NDBasis(6,i,k)+NDBasis(7,i,k).*x.*x+NDBasis(8,i,k).*x.*y;
curlBasis=@(x,y,i) NDBasis(2,i,k)-NDBasis(7,i,k).*x-2.*NDBasis(8,i,k).*y...
                 -NDBasis(4,i,k)-2.*NDBasis(7,i,k).*x-NDBasis(8,i,k).*y;

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
       f=@(x,y) basisX(x,y,s).*basisX(x,y,t)+basisY(x,y,s).*basisY(x,y,t)+curlBasis(x,y,s).*curlBasis(x,y,t);
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
   
    basisX=@(x,y,i) NDBasis(1,i,k).*x+NDBasis(2,i,k).*y+NDBasis(3,i,k)-NDBasis(7,i,k).*x.*y-NDBasis(8,i,k).*y.*y;
    basisY=@(x,y,i) NDBasis(4,i,k).*x+NDBasis(5,i,k).*y+NDBasis(6,i,k)+NDBasis(7,i,k).*x.*x+NDBasis(8,i,k).*x.*y;

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

fprintf("L2 Error mesh level %d: %e",mesh_level,error);
fprintf('\n');
