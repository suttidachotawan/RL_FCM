clear all
clc
filename = 'Unbalance.xlsx';
X = xlsread(filename);
er = 2.7182818;
r1=1;
r2=1;
r3=1;
t=1;
tol = 1e-9;
[nr nc] = size(X);
for i = 1:nc
   for j = 1:nr
     data(j, i) = (X(j, i)-std(X(:, i)))/mean(X(:, i));
    
   end
end
%data=X;
n = nr;
c_old = n;
c_new = n;
Max = n;
for i = 1:n
    A_old(i) = 1/n;
end
for i = 1:n
    
    A_new(i) = 1/n;
end
[n, no] = size(X);
c = n;
U = zeros([c, n]);
v = data;
U = rand([c, n]);
coutA_new = 0;
for i = 1:c
   for j = 1:n
     euclidean(j,i) = sqrt( ((data(j,1)-v(i,1)).^2) + ((data(j,2)-v(i,2)).^2) );
   end
end

v_old = v;
delta = 1e4;
k = 0;
while  (k<Max & delta>tol)
    
    for i = 1:nr
      for j = 1:n
        U(i, j)= exp(((-euclidean(j,i).^2)+(r1*log(A_new(i))))/r2)./sum(exp(((-euclidean(i,:).^2)+(r1*log(A_old(i))))/r2));
      end
    end
    r1 = er^(-1/10);
    r2 = er^(-1/100);
   
    for i =1:n
          SU(i) = sum(U(:,i)); 
          A_oldL(i) = log(A_new(i));
    end
    
        
  for i = 1:n
      A_new(i) = (1/n)*sum(SU(i))+((r3/r1)*A_new(i))*(log(A_new(i))-sum(A_oldL(i)));
  end
   % for i = 1:nr
     % for j = 1:no
       %checkn(i,j) = 2./(euclidean(i,:).^euclidean(i,:));
      %end
    %end
    N =0.004347053;
    %N=0.008347053;
    r3 = sum(exp((-N*n*(A_new(i)-A_new(i)))))/c;
    for i = 1:c_new
        if A_new(i) <= 1/n
          coutA_new = coutA_new+1;  
        end
    end 
    c_new = c_old - coutA_new;
    
    for j = 1:n
      U(:, j) = U(:, j)./sum(U(:, j)); 
    end  
 
     if t>=100 && c_old-c_new==0
        r3=0;
    else
      for i = 1:c
         v(i, :) = sum((data(:, :).*repmat(U(i, :)', 1, no)),1)./sum(U(i, :));
  
      end
    end
       
      for i = 1:c
          if v(i, :) > data(i,:)
              v(i, :) = 0;
              
          else
              if v(i, :)< data (i,:)
                  v(i, :) = 0;
              end
              v(i, :) = v(i, :);
          end
      end   
    


v_new = v;
delta = max(max(abs(v_new-v_old)));
v_old = v;
k = k+1;
t= t+1;
group = zeros([1, n]);

for i = 1:n
   
    [M, group(i)]=max(data(i, :));
end

end
     
for iter = 1:t
     iterations = iter
end

hold on
plot(data(1:6,1),data(1:6,2),'Ob')
plot(v(1:6,1),v(1:6,2),'xr','MarkerSize',10,'LineWidth',2)
hold off
