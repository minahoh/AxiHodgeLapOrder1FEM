function [RTBasis] = basis_functions_rt1(p,t,ele,edge,new_ele)
%BASIS_FUNCTIONS_RT1 Create a piecewise basis function for 
%   each edge and triangle of a triangulation
%
%    Structure of the variable "RTBasis". RTBasis(:,:,k) provides the information 
%    needed to get the eight basis functions for the k-th triangle. 
%    For i=1, 2, ..., 8, 
%    basisX=@(x,y,i) RTBasis(1,i,k).*x+RTBasis(2,i,k).*y+RTBasis(3,i,k)+RTBasis(7,i,k).*x.*x+RTBasis(8,i,k).*x.*y;
%    basisY=@(x,y,i) RTBasis(4,i,k).*x+RTBasis(5,i,k).*y+RTBasis(6,i,k)+RTBasis(7,i,k).*x.*y+RTBasis(8,i,k).*y.*y;
%    where basisX and basisY are the X and Y coordinates of the i-th basis
%    function respectively. 
%    Also not that,
%    divBasis=@(x,y,i) RTBasis(1,i,k)+2.*RTBasis(7,i,k).*x+RTBasis(8,i,k).*y...
%                  +RTBasis(5,i,k)+2.*RTBasis(8,i,k).*y+RTBasis(7,i,k).*x;
%
%    RTBasis(:,1:2,k) represents the two basis functions corresponding to the first edge of the k-th triangle.
%    RTBasis(:,3:4,k) represents the two basis functions corresponding to the second edge of the k-th triangle.
%    RTBasis(:,5:6,k) represents the two basis functions corresponding to the third edge of the k-th triangle.
%    RTBasis(:,7:8,k) represents the two basis functions corresponding to
%    the k-th triangle.
%
%    Global indexing of the basis functions: 
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
% Syntax:
%     basis = basis_functions_rt1(p,t,ele,edge,new_ele);
%
% Inputs:
%     p - a 2xNumNodes matrix representing nodal coordinates.
%     t - a 4xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs. The end row of T represents the geometry face ID 
%         to which the element belongs.
%     ele - a 3xNumTriangles matrix representing the element connectivity in 
%         terms of node IDs.
%     edge - a 2xNumEdges matrix representing each edge as a row with
%         starting node in column 1 and the ending node in column 2.
%     new_ele - a NumTrianglesx3 matrix representing the which edges
%         correspond to which triangles. t_e(T,i) represents the ith edge
%         in triangle T.
% Outputs:
%     basis - a matrix representing piece-wise basis functions for 
%         each edge and triangle in each triangle. basis_rt1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%
% Author: Minah Oh
% Edited By: Nicole Stock
% Date: Spring 2021


node=p';
[~,triangles]=size(t);

RTBasis=zeros(8,8,triangles);
abc=ones(8,8);

for k=1:triangles
   
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

      abc(2*j-1,1)=edge_length*(xe/2 + xs/2);
      abc(2*j-1,2)=edge_length*(ye/2 + ys/2);
      abc(2*j-1,3)=edge_length;
      abc(2*j-1,4)=0;
      abc(2*j-1,5)=0;
      abc(2*j-1,6)=0;
      abc(2*j-1,7)=edge_length*(xe^2/3 + (xe*xs)/3 + xs^2/3);
      abc(2*j-1,8)=edge_length*((xe*ye)/3 + (xe*ys)/6 + (xs*ye)/6 + (xs*ys)/3);

      abc(2*j,1)=edge_length*(xe/3 + xs/6);
      abc(2*j,2)=edge_length*(ye/3 + ys/6);
      abc(2*j,3)=edge_length/2;
      abc(2*j,4)=0;
      abc(2*j,5)=0;
      abc(2*j,6)=0;
      abc(2*j,7)=edge_length*(xe^2/4 + (xe*xs)/6 + xs^2/12);
      abc(2*j,8)=edge_length*((xe*ye)/4 + (xe*ys)/12 + (xs*ye)/12 + (xs*ys)/12);



    elseif ys==ye

     abc(2*j-1,1)=0;
     abc(2*j-1,2)=0;
     abc(2*j-1,3)=0;
     abc(2*j-1,4)=edge_length*(xe/2 + xs/2);
     abc(2*j-1,5)=edge_length*(ye/2 + ys/2);
     abc(2*j-1,6)=edge_length;
     abc(2*j-1,7)=edge_length*((xe*ye)/3 + (xe*ys)/6 + (xs*ye)/6 + (xs*ys)/3);
     abc(2*j-1,8)=edge_length*(ye^2/3 + (ye*ys)/3 + ys^2/3);

     abc(2*j,1)=0;
     abc(2*j,2)=0;
     abc(2*j,3)=0;
     abc(2*j,4)= edge_length*(xe/3 + xs/6);
     abc(2*j,5)=edge_length*(ye/3 + ys/6);
     abc(2*j,6)=edge_length/2;
     abc(2*j,7)=edge_length*((xe*ye)/4 + (xe*ys)/12 + (xs*ye)/12 + (xs*ys)/12);
     abc(2*j,8)=edge_length*(ye^2/4 + (ye*ys)/6 + ys^2/12);

    else

      abc(2*j-1,1)=edge_length*((2^(1/2)*xe)/4 + (2^(1/2)*xs)/4);
      abc(2*j-1,2)=edge_length*((2^(1/2)*ye)/4 + (2^(1/2)*ys)/4);
      abc(2*j-1,3)=(2^(1/2)*edge_length)/2;
      abc(2*j-1,4)=edge_length*((2^(1/2)*xe)/4 + (2^(1/2)*xs)/4);
      abc(2*j-1,5)=edge_length*((2^(1/2)*ye)/4 + (2^(1/2)*ys)/4);
      abc(2*j-1,6)=(2^(1/2)*edge_length)/2;
      abc(2*j-1,7)=edge_length*((2^(1/2)*xe^2)/6 + (2^(1/2)*xs^2)/6 + (2^(1/2)*xe*xs)/6 + (2^(1/2)*xe*ye)/6 + (2^(1/2)*xe*ys)/12 + (2^(1/2)*xs*ye)/12 + (2^(1/2)*xs*ys)/6);
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
 
end


end

