% clear all;
% close all;
% 
% load eigvals

N = length(eigvecorg)
ndim = length(eigvecorg(:,1,1))
gcs = 2


% ndim subplots, one for each of the eigenvalues.  Possibly fewer than 10 subplots needed because 
% there will usually be at least one or two complex pairs of eigenvalues.
hf=figure; 
rect = get(hf,'Position'); 
rect(1:2) = [0 0]; 
clear h

% Get initial axis handles and plot initial eigenvector components
for i = 1:4,
    subplot(2,2,i)
    hold on;
    h(1,i) = plot(real(eigvecorg(1,i,1)),imag(eigvecorg(1,i,1)),'ro');
    h(2,i) = plot([0;real(eigvecorg(1,i,1))],[0;imag(eigvecorg(1,i,1))],'r-')
    h(3,i) = plot(real(eigvecorg(2,i,1)),imag(eigvecorg(2,i,1)),'bx');
    h(4,i) = plot([0;real(eigvecorg(2,i,1))],[0;imag(eigvecorg(2,i,1))],'b-')
    s = sprintf('\\lambda = (%2.2f, %2.2f j)',real(eigorg(i,1)),imag(eigorg(i,1)))
    tit(i) = title(s)
%     legend('Lean','Steer',0)
    if i == 1,
        axis([-0.4,0.4,-0.4,0.4])
    elseif i == 2,
        axis([-0.4,0.4,-0.4,0.4])
    elseif i == 3,
        axis([-0.2,1,-1,1])
    elseif i == 4,
        axis([-0.2,0.05,-0.2,0.2])
%     axis([-0.9, 0.9,-0.9,0.9])
    end
    axis square
    box on;
    grid off
end



% Loop over each velocity
for k=2:N,
%     Loop over the first gcs entries of each vector, corresponding to each
%     eigenvalue
    for i = 1:4,
        subplot(2,2,i);
        XData = cell(2,1);
        YData = cell(2,1);
        for j=1:2,
            XData{j} = [0., real(eigvecorg(j,i,k))];
            YData{j} = [0., imag(eigvecorg(j,i,k))];
        end
        set(h(1,i),'XData',XData{1}(2))
        set(h(1,i),'YData',YData{1}(2))
        set(h(2,i),'XData',XData{1})
        set(h(2,i),'YData',YData{1})
        set(h(3,i),'XData',XData{2}(2))
        set(h(3,i),'YData',YData{2}(2))
        set(h(4,i),'XData',XData{2})
        set(h(4,i),'YData',YData{2})
        s = sprintf('\\lambda = (%2.2f, %2.2f j)',real(eigorg(i,k)),imag(eigorg(i,k)));
        set(tit(i),'String',s)
    end
    M(k) = getframe(hf,rect);
end


movie2avi(M(2:end),'eigvec.avi','compression','Indeo5','fps',15,'quality',100)
