function [Grad] = GradientEst(Node,Pot);

% Compute the gradient of the potential at each node 
% using a least square system.
% Where the potential of a Node is expressed at Taylor series approximation
% of each neighbor potentials
% f(x_1,y_1) = f(x_2,y_2) + grad f_x (x_2-x_1)+ grad f_x (y_2-y_1)
% INPUT
%
% Node = nodal structure
% Pot = potential values
% OUTPUT
%
% Grad = gradient Nx2 matrix



msN = max(size(Node));
dim = 2;%size(Node.Coordinate); 
Grad = zeros(msN,2);
%MaxEdgeL = 0;
for ii=1:msN
 Inds=[Node(ii).NodeConnection];
 NumCon = length(Inds);
 A = zeros(NumCon,dim);
 b = zeros(NumCon,1);
 W = ones(NumCon,NumCon);
 for i = 1: NumCon 
     
    dx = Node(ii).Coordinate(1,1) - Node(Inds(i)).Coordinate(1,1);
    dy = Node(ii).Coordinate(1,2) - Node(Inds(i)).Coordinate(1,2);
    A(i,:) = [dx dy];
    W(i,i) = 1/(dx*dx+dy*dy);
 
 end

 b = ones(NumCon,1)*Pot(ii) - Pot(Inds);
 
 Grad(ii,:) = ((A'*W^2*A)\(A'*W^2*b) )';  
end

%Plot gradient

% X=zeros(msN,1);
% Y=zeros(size(X));
% for i = 1:msN
%     X(i) = Node(i).Coordinate(1,1);
%     Y(i) = Node(i).Coordinate(1,2);
% end
%figure
%VectorFieldQuiver2D(X,Y,-Grad(:,1),-Grad(:,2),MaxEdgeL,10);

%hold on
%quiver(X,Y,Grad(:,1),Grad(:,2),'r')

%hold off