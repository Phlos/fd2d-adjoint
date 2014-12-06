% First, create 4  figures with four different graphs (each with  a 

% colorbar):

figure(1)

surf(peaks(10))

colorbar

figure(2)

mesh(peaks(10))

colorbar

figure(3)

contour(peaks(10))

colorbar

figure(4)

pcolor(peaks(10))

colorbar

% Now create destination graph

figure(5)

ax = zeros(4,1);

% for i = 1:4
% 
%     ax(i)=subplot(4,1,i);
% 
% end

ax(1) = subplot(4,4,[2,4]-1);
ax(2) = subplot(4,4,3+[2,4]);
ax(3) = subplot(4,4,7+[2,4]);
ax(4) = subplot(4,4,11+[2,4]);

pause(5);

% Now copy contents of each figure over to destination figure

% Modify position of each axes as it is transferred

for i = 1:4

    figure(i)

    h = get(gcf,'Children');

    newh = copyobj(h,5);

    for j = 1:length(newh)

posnewh = get(newh(j),'Position');

possub  = get(ax(i),'Position');

set(newh(j),'Position', ...
[posnewh(1) possub(2) posnewh(3) possub(4)])

    end

    delete(ax(i));

end

figure(5)
