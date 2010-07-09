clear all
close all
load bfr_eig_data
eigvalreal = real(eigval); % extract real parts of eigenvalues
eigvalimag = imag(eigval); % extract the imaginary parts
figure(1)
hold on
N=length(v)
% plot(v,eigvalreal( 1,1:N),'.b')
% plot(v,eigvalreal( 2,1:N),'og')
% plot(v,eigvalreal( 3,1:N),'.r')
% plot(v,eigvalreal( 4,1:N),'oc')
% plot(v,eigvalreal( 5,1:N),'.m')
% plot(v,eigvalreal( 6,1:N),'oy')
% plot(v,eigvalreal( 7,1:N),'.k')
% plot(v,eigvalreal( 8,1:N),'ob')
% plot(v,eigvalreal( 9,1:N),'.g')
% plot(v,eigvalreal(10,1:N),'or')
% plot(v,zeros(N,1),'k')
% hold off
%-------------------Organize Eigenvalue Matrix
% set first column in organized matrix to the first column of the original matrix
eigorg(:,1)=eigval(:,1);
% rearrange the eigenvalue columns by calcuating the absolute value of the
% difference between values of the successicve column and the preceeding
% column

for i = 2
    dist=zeros(10,10);
    for j=1:10,
        dist(i,j)=abs(eigval(i,j)-eigval(i-1,j));
    end
end
%-------------------Plot the Eigenvalues
% eigvalreal = real(eigorg); % extract real parts of eigenvalues
% eigvalimag = imag(eigorg); % extract the imaginary parts
% figure(2)
% hold on
% plot(v,eigvalreal( 1,1:N),'.b')
% plot(v,eigvalreal( 2,1:N),'og')
% plot(v,eigvalreal( 3,1:N),'.r')
% plot(v,eigvalreal( 4,1:N),'oc')
% plot(v,eigvalreal( 5,1:N),'.m')
% plot(v,eigvalreal( 6,1:N),'oy')
% plot(v,eigvalreal( 7,1:N),'.k')
% plot(v,eigvalreal( 8,1:N),'ob')
% plot(v,eigvalreal( 9,1:N),'.g')
% plot(v,eigvalreal(10,1:N),'or')
% plot(v,zeros(N,1),'k')
% hold off
% 
% 
% 






% 
% for i = 1:length(v)-1
%     indiceused=zeros(10,1);
%     for j = 1:10
%         first = eigorg(j,i);
%         for k =1:10
%             second = eigval(k,i+1);
%             diff(k) = abs(second - first);
%         end
%         [mindiff indice]=min(diff);
%         for m=1:length(indiceused)
%             x=indice-indiceused(m);
%             if x == 0
%                 check=1; % the indice has been used!
%             end
%         end
%         if check ~= 1 % if the indice hasn't been used
%             eigorg(j,i+1)=eigval(indice,i+1);
%         else % if the indice has been used
%             maxdiff=max(diff);
%             diff(indice)=max(diff)+1;
%             [mindiff indice]=min(diff);
%             eigorg(j,i+1)=eigval(indice,i+1);
%         end
%         indiceused(j)=indice;
%     end
% end
% 



