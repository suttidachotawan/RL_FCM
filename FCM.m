clear all
clc
filename = 'Unbalance.xlsx';
X = xlsread(filename);
m = 1.5;
Max = 1000;
tol = 1e-3;
c = 8;

[nr nc] = size(X);
for i = 1:nc
   for j = 1:nr
     data(j, i) = (X(j, i)-std(X(:, i)))/mean(X(:, i));
    
   end
end

[n, no] = size(X);
U = zeros([c, n]);
v = repmat(max(data), c, 1).*rand([c, no]);
U = rand([c, n]);
randomcenter = v

for i = 1:c
   for j = 1:n
     euclidean(j,i) = sqrt( (data(j,1)-v(i,1)).^2 + (data(j,2)-v(i,2)).^2 );
   end
end

for j = 1:n
      U(:, j) = U(:, j)./sum(U(:, j));         
end  

for i = 1:c
      v(i, :) = sum((data(:, :).*repmat(U(i, :).^m', 1, no)),1)./sum(U(i, :).^m);
end
v_old = v;
delta = 1e4;
k = 0;
Color=0;
while  (k<Max & delta>tol)
    
    for i = 1:c
      for j = 1:n
        U(i, j) = (1./euclidean(j,i)).^(1/(m-1))/sum(1./euclidean(j,:).^(1/(m-1)));
       
      end
    end
    for i = 1:c
        v(i,:) = sum((data(:, :).*repmat(U(i, :).^m', 1, no)), 1)./sum(U(i, :).^m);
    end

for i = 1:c
   for j = 1:n
     euclidean(j,i) = sqrt( (data(j,1)-v(i,1)).^2 + (data(j,2)-v(i,2)).^2 );
   end
end
v_new = v;
delta = max(max(abs(v_new-v_old)));
v_old = v;
k = k+1;

prediction = zeros([n, 1]);
for i = 1:n
   [M, prediction(i)]=max(U(:, i));
  
end
end
iterations1 = k;

for iter = 1:iterations1
     iterations = iter
end

hold on
xlabel('Price USD');
ylabel('Marketcap USD'); 
plot(data(1:1,1),data(1:1,2),'ok');
plot(data(2:2,1),data(2:2,2),'ob');
plot(data(3:4,1),data(3:4,2),'og');
plot(data(5:500,1),data(5:500,2),'om');
%plot(data(1:240,1),data(1:240,2),'Ob')
%plot(v(1:c,1),v(1:c,2),'xr','MarkerSize',10,'LineWidth',2)
hold off
