function plotAbundance(abudance,algorithm,cood)
figure;
massNum=6;
width=64;
len=64;
for i=1:massNum
    temp=abudance(:,i);
    temp=temp-min(temp);
    temp=uint8(temp/max(temp)*255);
    temp=reshape(temp,[width len]);
%     subplot(1,massNum,i);
    imshow(temp);
    saveas(gcf,strcat(algorithm,'Abundance',cood{i},'.eps'));
    close gcf;
%     title('Estimated');
%     subplot(2,massNum,i+massNum);
%     temp=trueAbundance(:,i);
%     temp=temp-min(temp);
%     temp=uint8(temp/max(temp)*255);
%     temp=reshape(temp,[width len]);
%     imshow(temp);
%     title('Ground-Truth');
end