function [ClippedImage]=ClipImage(InputImage)
[m n]=size(InputImage);
min_row=0;
min_col=0;
for i=1:m
  if(sum(InputImage(i,:))==0 && min_row~=0)
      max_row=i-1;
      break;
  end
  if(sum(InputImage(i,:)) >= 1 && min_row == 0)
      min_row = i;
  end
end

for i=1:n
   if(sum(InputImage(:,i))==0 && min_col~=0)
      max_col=i-1;
      break;
  end
  if(sum(InputImage(:,i)) >= 1 && min_col == 0)
      min_col = i;
  end
end 

ClippedImage=InputImage(min_row:max_row,min_col:max_col);
