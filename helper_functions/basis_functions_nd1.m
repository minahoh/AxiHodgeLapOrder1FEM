function [NDBasis] = basis_functions_nd1(p,t,ele,edge,new_ele)
%BASIS_FUNCTIONS_RT1 Create a piecewise basis function for 
%   each edge and triangle of a triangulation
%
%    Structure of the variable "NDBasis". NDBasis(:,:,k) provides the information 
%    needed to get the eight basis functions for the k-th triangle. 
%    For i=1, 2, ..., 8, 
%    basisX=@(x,y,i) NDBasis(1,i,k).*x+NDBasis(2,i,k).*y+NDBasis(3,i,k)-NDBasis(7,i,k).*x.*y-NDBasis(8,i,k).*y.*y;
%    basisY=@(x,y,i) NDBasis(4,i,k).*x+NDBasis(5,i,k).*y+NDBasis(6,i,k)+NDBasis(7,i,k).*x.*x+NDBasis(8,i,k).*x.*y;
%    where basisX and basisY are the X and Y coordinates of the i-th basis
%    function respectively. 
%
%    NDBasis(:,1:2,k) represents the two basis functions corresponding to the first edge of the k-th triangle.
%    NDBasis(:,3:4,k) represents the two basis functions corresponding to the second edge of the k-th triangle.
%    NDBasis(:,5:6,k) represents the two basis functions corresponding to the third edge of the k-th triangle.
%    NDBasis(:,7:8,k) represents the two basis functions corresponding to
%    the k-th triangle.
%
%    Global indexing of the basis functions: 
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
%         each edge and triangle in each triangle. basis_nd1(:,i,T) represents the 
%         pieceiwise basis function for the ith edge in triangle T.
%
% Author: Minah Oh
% Edited By: Nicole Stock
% Date: Spring 2021

node=p';
[~,triangles]=size(t);

NDBasis=zeros(8,8,triangles);
abc=ones(8,8);

for k=1:triangles

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
end