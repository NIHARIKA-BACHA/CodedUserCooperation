function [G Q1 Q2 j_limit k_limit] = ldgm_G(M,j,k)
%   M= 12;
%   j= 3;
%   k= 8;

P=zeros(k,M);  %generate the empty P matrix  
  
%number of parity bits
%parity_bits = N*j/k;
%group = M/k;
%==========================================================================
                %Generating parity check matrix H
 P = generate(); % Generate P matrix
 
 function P = generate()
   for i=1:M % Making of parity matrix
    column = [ones(1,j) zeros(1,k-j)];
    column = column(randperm(k)); % Random permutation of edges (1s)
    P(:,i) = column'; % Final P matrix
   end
   %P
   for i=1:k
    if length(find(P(i,:))) < 2 % Check that every bit node is connected to at least 2 check nodes
       P = generate();
    end
   end
 end
%P
%G = cat(2, eye(k),P);
G = P;

H = [P' eye(M)];

Q1=[];
Q2=[];
Q1row = 1;
for col_index = 1:M+k
    
    for row_index = 1:M
        if H(row_index, col_index) == 1
        Q1(Q1row,col_index)= row_index;
        Q1row = Q1row + 1;
        end
    end
    Q1row = 1;
end
%Q1
Q2col=1;
Q2row = 1;
for row_index = 1:M
    
    for col_index = 1:M+k
        
        if H(row_index, col_index) == 1
        Q2(Q2row,Q2col)= col_index;
        Q2col = Q2col + 1;
        end
    end
    Q2col = 1;
    Q2row = Q2row +1;
end
%Q2

%================= Computing the j (1s in each column) and k (1s in each
%row) for LT Generator matrix ==========================================
for x1=1:M+k
    j_limit(x1)=length(find(nonzeros(Q1(:,x1))));
end
for x1=1:M
    k_limit(x1)=length(find(nonzeros(Q2(x1,:))));
end
%j_limit
%k_limit
end